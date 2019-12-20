#  ___________________________________________________________________________
#
#  Pyomo: Python Optimization Modeling Objects
#  Copyright 2017 National Technology and Engineering Solutions of Sandia, LLC
#  Under the terms of Contract DE-NA0003525 with National Technology and
#  Engineering Solutions of Sandia, LLC, the U.S. Government retains certain
#  rights in this software.
#  This software is distributed under the 3-clause BSD License.
#  ___________________________________________________________________________

#
# Problem Writer for CPLEX LP Format Files
#

import logging
import math
import operator

from six import iterkeys, iteritems, StringIO
from six.moves import xrange

from pyutilib.misc import PauseGC

from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.opt import ProblemFormat
from pyomo.opt.base import AbstractProblemWriter, WriterFactory
from pyomo.core.base import \
    (SymbolMap, TextLabeler,
     NumericLabeler, Constraint, SortComponents,
     Var, value,
     SOSConstraint, Objective,
     ComponentMap, is_fixed)
from pyomo.repn import generate_standard_repn

logger = logging.getLogger('pyomo.core')

def _no_negative_zero(val):
    """Make sure -0 is never output. Makes diff tests easier."""
    if val == 0:
        return 0
    return val

def _get_bound(exp):
    if exp is None:
        return None
    if is_fixed(exp):
        return value(exp)
    raise ValueError("non-fixed bound or weight: " + str(exp))


@WriterFactory.register('new_cpxlp', 'Generate the corresponding CPLEX LP file')
@WriterFactory.register('new_lp', 'Generate the corresponding CPLEX LP file')
class ProblemWriter_cpxlp(AbstractProblemWriter):

    CONFIG = ConfigBlock()
    CONFIG.declare('symbolic_solver_labels', ConfigValue(
        default=False,
        domain=bool,
        doc="""Use full Pyomo component names in the LP file rather than
        shortened symbols (slower, but useful for debugging)."""
    ))
    CONFIG.declare('labeler', ConfigValue(
        default=None,
    ))
    CONFIG.declare('skip_trivial_constraints', ConfigValue(
        default=False,
        domain=bool,
        doc="""Skip writing constraints whose body section is fixed (i.e., no
        variables)"""
    ))
    CONFIG.declare('output_fixed_variable_bounds', ConfigValue(
        default=False,
        domain=bool,
    ))
    CONFIG.declare('include_all_variable_bounds', ConfigValue(
        default=False,
        domain=bool,
        doc="""If False, unused variables will not be included in the
        LP file. Otherwise, include all variables in the bounds sections."""
    ))
    CONFIG.declare('file_determinism', ConfigValue(
        default=1,
        domain=In([0,1,2]),
        doc="""Specify how much effort do we want to put into ensuring the
        LP file is written deterministically for a Pyomo model:
          - 0 : None
          - 1 : sort keys of indexed components (default)
          - 2 : sort keys AND sort names (over declaration order)
        """
    ))
    CONFIG.declare('row_order', ConfigValue(
        default=None,
        domain=list,
        doc="Specify the order that constraints appear in the output"
    ))
    CONFIG.declare('column_order', ConfigValue(
        default=None,
        domain=list,
        doc="Specify the order that variables appear in the output"
    ))
    CONFIG.declare('force_objective_constant', ConfigValue(
        default=False,
        domain=bool,
        doc="""Make sure the ONE_VAR_CONSTANT variable appears in the
        objective even if the constant part of the objective is zero"""
    ))
    CONFIG.declare('allow_quadratic_objective', ConfigValue(
        default=True,
        domain=bool,
    ))
    CONFIG.declare('allow_quadratic_constraints', ConfigValue(
        default=True,
        domain=bool,
    ))
    CONFIG.declare('allow_sos1', ConfigValue(
        default=True,
        domain=bool,
    ))
    CONFIG.declare('allow_sos2', ConfigValue(
        default=True,
        domain=bool,
    ))

    def __init__(self):

        AbstractProblemWriter.__init__(self, ProblemFormat.cpxlp)

        # The LP writer tracks which variables are
        # referenced in constraints, so that a user does not end up with a
        # zillion "unreferenced variables" warning messages.
        # This dictionary maps id(_VarData) -> _VarData.
        self._referenced_variable_ids = {}

        # Per ticket #4319, we are using %.17g, which mocks the
        # behavior of using %r (i.e., float('%r'%<number>) == <number>) with
        # the added benefit of outputting (+/-). The only case where this
        # fails to mock the behavior of %r is for large (long) integers (L),
        # which is a rare case to run into and is probably indicative of
        # other issues with the model.
        # *** NOTE ***: If you use 'r' or 's' here, it will break code that
        #               relies on using '%+' before the formatting character
        #               and you will need to go add extra logic to output
        #               the number's sign.
        self._precision_string = '.17g'
        self.linear_coef_string_template = '%+'+self._precision_string+' %s\n'
        self.quad_coef_string_template =   '%+'+self._precision_string+' '
        self.obj_string_template =         '%+'+self._precision_string+' %s\n'
        self.sos_template_string =         "%s:%"+self._precision_string+"\n"
        self.eq_string_template =          "= %"+self._precision_string+'\n'
        self.geq_string_template =         ">= %"+self._precision_string+'\n\n'
        self.leq_string_template =         "<= %"+self._precision_string+'\n\n'
        self.lb_string_template =          "%"+self._precision_string+" <= "
        self.ub_string_template =          " <= %"+self._precision_string+"\n"

    def __call__(self,
                 model,
                 output_filename,
                 io_options=None,
                 **kwds):

        config = self.CONFIG(io_options)(kwds)

        if config.symbolic_solver_labels and (config.labeler is not None):
            raise ValueError("ProblemWriter_cpxlp: Using both the "
                             "'symbolic_solver_labels' and 'labeler' "
                             "I/O options is forbidden")

        #
        # Create labeler
        #
        if config.symbolic_solver_labels:
            labeler = TextLabeler()
        elif config.labeler is not None:
            labeler = config.labeler
        else:
            labeler = NumericLabeler('x')

        if output_filename is None:
            output_filename = model.name + ".lp"

        # when sorting, there are a non-trivial number of temporary
        # objects created. these all yield non-circular references, so
        # disable GC - the overhead is non-trivial, and because
        # references are non-circular, in CPython everything will be
        # collected immediately anyway.
        with PauseGC() as pgc:
            with open(output_filename, "w") as output_file:
                symbol_map = self._print_model_LP(
                    config, model, output_file, labeler)

        return output_filename, symbol_map

    def _print_expr_canonical(self,
                              x,
                              output,
                              object_symbol_dictionary,
                              variable_symbol_dictionary,
                              is_objective,
                              config):

        """
        Return a expression as a string in LP format.

        Note that this function does not handle any differences in LP format
        interpretation by the solvers (e.g. CPlex vs GLPK).  That decision is
        left up to the caller.

        required arguments:
          x: A Pyomo canonical expression to write in LP format
        """
        assert (not config.force_objective_constant) or (is_objective)
        linear_coef_string_template = self.linear_coef_string_template
        quad_coef_string_template = self.quad_coef_string_template

        constant=True
        #
        # Linear
        #
        if len(x.linear_vars) > 0:
            constant=False
            for vardata in x.linear_vars:
                self._referenced_variable_ids[id(vardata)] = vardata

            if config.column_order is None:
                #
                # Order columns by dictionary names
                #
                names = [variable_symbol_dictionary[id(var)] for var in x.linear_vars]

                for i, name in sorted(enumerate(names), key=lambda x: x[1]):
                    output.append(linear_coef_string_template % (x.linear_coefs[i], name))
            else:
                #
                # Order columns by the value of column_order[]
                #
                for i, var in sorted(enumerate(x.linear_vars), key=lambda x: config.column_order[x[1]]):
                    name = variable_symbol_dictionary[id(var)]
                    output.append(linear_coef_string_template % (x.linear_coefs[i], name))
        #
        # Quadratic
        #
        if len(x.quadratic_vars) > 0:
            constant=False
            for var1, var2 in x.quadratic_vars:
                self._referenced_variable_ids[id(var1)] = var1
                self._referenced_variable_ids[id(var2)] = var2

            output.append("+ [\n")

            if config.column_order is None:
                #
                # Order columns by dictionary names
                #
                quad = set()
                names = []
                i = 0
                for var1, var2 in x.quadratic_vars:
                    name1 = variable_symbol_dictionary[id(var1)]
                    name2 = variable_symbol_dictionary[id(var2)]
                    if name1 < name2:
                        names.append( (name1,name2) )
                    elif name1 > name2:
                        names.append( (name2,name1) )
                    else:
                        quad.add(i)
                        names.append( (name1,name1) )
                    i += 1
                for i, names_ in sorted(enumerate(names), key=lambda x: x[1]):
                    #
                    # Times 2 because LP format requires /2 for all the quadratic
                    # terms /of the objective only/.  Discovered the last bit thru
                    # trial and error.  Obnoxious.
                    # Ref: ILog CPlex 8.0 User's Manual, p197.
                    #
                    if is_objective:
                        tmp = 2*x.quadratic_coefs[i]
                        output.append(quad_coef_string_template % tmp)
                    else:
                        output.append(quad_coef_string_template % x.quadratic_coefs[i])
                    if i in quad:
                        output.append("%s ^ 2\n" % (names_[0]))
                    else:
                        output.append("%s * %s\n" % (names_[0], names_[1]))
            else:
                #
                # Order columns by the value of column_order[]
                #
                quad = set()
                cols = []
                i = 0
                for var1, var2 in x.quadratic_vars:
                    col1 = config.column_order[var1]
                    col2 = config.column_order[var2]
                    if col1 < col2:
                        cols.append( (((col1,col2) , variable_symbol_dictionary[id(var1)], variable_symbol_dictionary[id(var2)])) )
                    elif col1 > col2:
                        cols.append( (((col2,col1) , variable_symbol_dictionary[id(var2)], variable_symbol_dictionary[id(var1)])) )
                    else:
                        quad.add(i)
                        cols.append( ((col1,col1), variable_symbol_dictionary[id(var1)]) )
                    i += 1
                for i, cols_ in sorted(enumerate(cols), key=lambda x: x[1][0]):
                    #
                    # Times 2 because LP format requires /2 for all the quadratic
                    # terms /of the objective only/.  Discovered the last bit thru
                    # trial and error.  Obnoxious.
                    # Ref: ILog CPlex 8.0 User's Manual, p197.
                    #
                    if is_objective:
                        output.append(quad_coef_string_template % 2*x.quadratic_coefs[i])
                    else:
                        output.append(quad_coef_string_template % x.quadratic_coefs[i])
                    if i in quad:
                        output.append("%s ^ 2\n" % cols_[1])
                    else:
                        output.append("%s * %s\n" % (cols_[1], cols_[2]))

            output.append("]")

            if is_objective:
                output.append(' / 2\n')
                # divide by 2 because LP format requires /2 for all the quadratic
                # terms.  Weird.  Ref: ILog CPlex 8.0 User's Manual, p197
            else:
                output.append("\n")

        if constant and not is_objective:
            # If we made it to here we are outputing
            # trivial constraints place 0 *
            # ONE_VAR_CONSTANT on this side of the
            # constraint for the benefit of solvers like
            # Glpk that cannot parse an LP file without
            # a variable on the left hand side.
            output.append(linear_coef_string_template % (0, 'ONE_VAR_CONSTANT'))

        #
        # Constant offset
        #
        # Currently, it appears that we only need to print
        # the constant offset term for objectives.
        #
        if is_objective and (config.force_objective_constant or (x.constant != 0.0)):
            output.append(self.obj_string_template % (x.constant, 'ONE_VAR_CONSTANT'))

        #
        # Return constant offset
        #
        return x.constant

    def printSOS(self,
                 symbol_map,
                 labeler,
                 variable_symbol_map,
                 soscondata,
                 output):
        """
        Prints the SOS constraint associated with the _SOSConstraintData object
        """
        sos_template_string = self.sos_template_string

        if hasattr(soscondata, 'get_items'):
            sos_items = list(soscondata.get_items())
        else:
            sos_items = list(soscondata.items())

        if len(sos_items) == 0:
            return

        level = soscondata.level

        output.append('%s: S%s::\n'
                          % (symbol_map.getSymbol(soscondata,labeler), level))

        for vardata, weight in sos_items:
            weight = _get_bound(weight)
            if weight < 0:
                raise ValueError(
                    "Cannot use negative weight %f "
                    "for variable %s is special ordered "
                    "set %s " % (weight, vardata.name, soscondata.name))
            if vardata.fixed:
                raise RuntimeError(
                    "SOSConstraint '%s' includes a fixed variable '%s'. This is "
                    "currently not supported. Deactive this constraint in order to "
                    "proceed." % (soscondata.name, vardata.name))
            self._referenced_variable_ids[id(vardata)] = vardata
            output.append(sos_template_string
                              % (variable_symbol_map.getSymbol(vardata),
                                 weight))

    def _print_model_LP(self,
                        config,
                        model,
                        output_file,
                        labeler):

        eq_string_template = self.eq_string_template
        leq_string_template = self.leq_string_template
        geq_string_template = self.geq_string_template
        ub_string_template = self.ub_string_template
        lb_string_template = self.lb_string_template

        symbol_map = SymbolMap()
        variable_symbol_map = SymbolMap()
        # NOTE: we use createSymbol instead of getSymbol because we
        #       know whether or not the symbol exists, and don't want
        #       to the overhead of error/duplicate checking.
        # cache frequently called functions
        create_symbol_func = SymbolMap.createSymbol
        create_symbols_func = SymbolMap.createSymbols
        alias_symbol_func = SymbolMap.alias
        variable_label_pairs = []

        # populate the symbol map in a single pass.
        #objective_list, constraint_list, sosconstraint_list, variable_list \
        #    = self._populate_symbol_map(model,
        #                                symbol_map,
        #                                labeler,
        #                                variable_symbol_map,
        #                                file_determinism=file_determinism)
        sortOrder = SortComponents.unsorted
        if config.file_determinism >= 1:
            sortOrder = sortOrder | SortComponents.indices
            if config.file_determinism >= 2:
                sortOrder = sortOrder | SortComponents.alphabetical

        #
        # Create variable symbols (and cache the block list)
        #
        all_blocks = []
        variable_list = []
        #
        # WEH - TODO:  See if this is faster
        # NOTE: This loop doesn't find all of the variables.  :(
        #
        #for block in model.block_data_objects(active=True,
        #                                      sort=sortOrder):
        #
        #    all_blocks.append(block)
        #
        #    for vardata in block.component_data_objects(
        #            Var,
        #            active=True,
        #            sort=sortOrder,
        #            descend_into=False):
        #
        #        variable_list.append(vardata)
        #        variable_label_pairs.append(
        #            (vardata,create_symbol_func(symbol_map,
        #                                        vardata,
        #                                        labeler)))
        all_blocks = list( model.block_data_objects(
                active=True, sort=sortOrder) )
        variable_list = list( model.component_data_objects(
                Var, sort=sortOrder) )
        variable_label_pairs = list(
            (vardata, create_symbol_func(symbol_map, vardata, labeler))
            for vardata in variable_list )
        variable_symbol_map.addSymbols(variable_label_pairs)

        # and extract the information we'll need for rapid labeling.
        object_symbol_dictionary = symbol_map.byObject
        variable_symbol_dictionary = variable_symbol_map.byObject

        # cache - these are called all the time.
        print_expr_canonical = self._print_expr_canonical

        # print the model name and the source, so we know roughly where
        # it came from.
        #
        # NOTE: this *must* use the "\* ... *\" comment format: the GLPK
        # LP parser does not correctly handle other formats (notably, "%").
        output = []
        output.append(
            "\\* Source Pyomo model name=%s *\\\n\n" % (model.name,) )

        #
        # Objective
        #

        supports_quadratic_objective = config.allow_quadratic_objective

        numObj = 0
        onames = []
        for block in all_blocks:

            gen_obj_repn = getattr(block, "_gen_obj_repn", True)

            # Get/Create the ComponentMap for the repn
            if not hasattr(block,'_repn'):
                block._repn = ComponentMap()
            block_repn = block._repn

            for objective_data in block.component_data_objects(
                    Objective,
                    active=True,
                    sort=sortOrder,
                    descend_into=False):

                numObj += 1
                onames.append(objective_data.name)
                if numObj > 1:
                    raise ValueError(
                        "More than one active objective defined for input "
                        "model '%s'; Cannot write legal LP file\n"
                        "Objectives: %s" % (model.name, ' '.join(onames)))

                create_symbol_func(symbol_map,
                                   objective_data,
                                   labeler)

                symbol_map.alias(objective_data, '__default_objective__')
                if objective_data.is_minimizing():
                    output.append("min \n")
                else:
                    output.append("max \n")

                if gen_obj_repn:
                    repn = generate_standard_repn(objective_data.expr)
                    block_repn[objective_data] = repn
                else:
                    repn = block_repn[objective_data]

                degree = repn.polynomial_degree()

                if degree == 0:
                    logger.warning("Constant objective detected, replacing "
                          "with a placeholder to prevent solver failure.")
                    config.force_objective_constant = True
                elif degree == 2:
                    if not supports_quadratic_objective:
                        raise RuntimeError(
                            "Selected solver is unable to handle "
                            "objective functions with quadratic terms. "
                            "Objective at issue: %s."
                            % objective_data.name)
                elif degree is None:
                    raise RuntimeError(
                        "Cannot write legal LP file.  Objective '%s' "
                        "has nonlinear terms that are not quadratic."
                        % objective_data.name)

                output.append(
                    object_symbol_dictionary[id(objective_data)]+':\n')

                offset = print_expr_canonical(
                    repn,
                    output,
                    object_symbol_dictionary,
                    variable_symbol_dictionary,
                    True,
                    config)

        if numObj == 0:
            raise ValueError(
                "ERROR: No objectives defined for input model. "
                "Cannot write legal LP file.")

        # Constraints
        #
        # If there are no non-trivial constraints, you'll end up with an empty
        # constraint block. CPLEX is OK with this, but GLPK isn't. And
        # eliminating the constraint block (i.e., the "s.t." line) causes GLPK
        # to whine elsewhere. Output a warning if the constraint block is empty,
        # so users can quickly determine the cause of the solve failure.

        output.append("\n")
        output.append("s.t.\n")
        output.append("\n")

        have_nontrivial = False

        supports_quadratic_constraint = config.allow_quadratic_constraints

        def constraint_generator():
            for block in all_blocks:

                gen_con_repn = getattr(block, "_gen_con_repn", True)

                # Get/Create the ComponentMap for the repn
                if not hasattr(block,'_repn'):
                    block._repn = ComponentMap()
                block_repn = block._repn

                for constraint_data in block.component_data_objects(
                        Constraint,
                        active=True,
                        sort=sortOrder,
                        descend_into=False):

                    if (not constraint_data.has_lb()) and \
                       (not constraint_data.has_ub()):
                        assert not constraint_data.equality
                        continue # non-binding, so skip

                    if constraint_data._linear_canonical_form:
                        repn = constraint_data.canonical_form()
                    elif gen_con_repn:
                        repn = generate_standard_repn(constraint_data.body)
                        block_repn[constraint_data] = repn
                    else:
                        repn = block_repn[constraint_data]

                    yield constraint_data, repn

        if config.row_order is not None:
            sorted_constraint_list = list(constraint_generator())
            sorted_constraint_list.sort(key=lambda x: config.row_order[x[0]])
            def yield_all_constraints():
                for data, repn in sorted_constraint_list:
                    yield data, repn
        else:
            yield_all_constraints = constraint_generator

        # FIXME: This is a hack to get nested blocks working...
        for constraint_data, repn in yield_all_constraints():
            have_nontrivial = True

            degree = repn.polynomial_degree()

            #
            # Write constraint
            #

            # There are conditions, e.g., when fixing variables, under which
            # a constraint block might be empty.  Ignore these, for both
            # practical reasons and the fact that the CPLEX LP format
            # requires a variable in the constraint body.  It is also
            # possible that the body of the constraint consists of only a
            # constant, in which case the "variable" of
            if degree == 0:
                if config.skip_trivial_constraints:
                    continue
            elif degree == 2:
                if not supports_quadratic_constraint:
                    raise ValueError(
                        "Solver unable to handle quadratic expressions. Constraint"
                        " at issue: '%s'" % (constraint_data.name))
            elif degree is None:
                raise ValueError(
                    "Cannot write legal LP file.  Constraint '%s' has a body "
                    "with nonlinear terms." % (constraint_data.name))

            # Create symbol
            con_symbol = create_symbol_func(symbol_map, constraint_data, labeler)

            if constraint_data.equality:
                assert value(constraint_data.lower) == \
                    value(constraint_data.upper)
                label = 'c_e_%s_' % con_symbol
                alias_symbol_func(symbol_map, constraint_data, label)
                output.append(label)
                output.append(':\n')
                offset = print_expr_canonical(repn,
                                              output,
                                              object_symbol_dictionary,
                                              variable_symbol_dictionary,
                                              False,
                                              config)
                bound = constraint_data.lower
                bound = _get_bound(bound) - offset
                output.append(eq_string_template
                                  % (_no_negative_zero(bound)))
                output.append("\n")
            else:
                if constraint_data.has_lb():
                    if constraint_data.has_ub():
                        label = 'r_l_%s_' % con_symbol
                    else:
                        label = 'c_l_%s_' % con_symbol
                    alias_symbol_func(symbol_map, constraint_data, label)
                    output.append(label)
                    output.append(':\n')
                    offset = print_expr_canonical(repn,
                                                  output,
                                                  object_symbol_dictionary,
                                                  variable_symbol_dictionary,
                                                  False,
                                                  config)
                    bound = constraint_data.lower
                    bound = _get_bound(bound) - offset
                    output.append(geq_string_template
                                      % (_no_negative_zero(bound)))
                else:
                    assert constraint_data.has_ub()

                if constraint_data.has_ub():
                    if constraint_data.has_lb():
                        label = 'r_u_%s_' % con_symbol
                    else:
                        label = 'c_u_%s_' % con_symbol
                    alias_symbol_func(symbol_map, constraint_data, label)
                    output.append(label)
                    output.append(':\n')
                    offset = print_expr_canonical(repn,
                                                  output,
                                                  object_symbol_dictionary,
                                                  variable_symbol_dictionary,
                                                  False,
                                                  config)
                    bound = constraint_data.upper
                    bound = _get_bound(bound) - offset
                    output.append(leq_string_template
                                      % (_no_negative_zero(bound)))
                else:
                    assert constraint_data.has_lb()

            # A simple hack to avoid caching super large files
            if len(output) > 1024:
                output_file.write( "".join(output) )
                output = []

        if not have_nontrivial:
            logger.warning('Empty constraint block written in LP format '  \
                  '- solver may error')

        # the CPLEX LP format doesn't allow constants in the objective (or
        # constraint body), which is a bit silly.  To avoid painful
        # book-keeping, we introduce the following "variable", constrained
        # to the value 1.  This is used when quadratic terms are present.
        # worst-case, if not used, is that CPLEX easily pre-processes it out.
        prefix = ""
        output.append('%sc_e_ONE_VAR_CONSTANT: \n' % prefix)
        output.append('%sONE_VAR_CONSTANT = 1.0\n' % prefix)
        output.append("\n")

        # SOS constraints
        #
        # For now, we write out SOS1 and SOS2 constraints in the cplex format
        #
        # All Component objects are stored in model._component, which is a
        # dictionary of {class: {objName: object}}.
        #
        # Consider the variable X,
        #
        #   model.X = Var(...)
        #
        # We print X to CPLEX format as X(i,j,k,...) where i, j, k, ... are the
        # indices of X.
        #
        SOSlines = []
        sos1 = config.allow_sos1
        sos2 = config.allow_sos2
        writtenSOS = False
        for block in all_blocks:

            for soscondata in block.component_data_objects(
                    SOSConstraint,
                    active=True,
                    sort=sortOrder,
                    descend_into=False):

                create_symbol_func(symbol_map, soscondata, labeler)

                level = soscondata.level
                if (level == 1 and not sos1) or \
                   (level == 2 and not sos2) or \
                   (level > 2):
                    raise ValueError(
                        "Solver does not support SOS level %s constraints" % (level))
                if writtenSOS == False:
                    SOSlines.append("SOS\n")
                    writtenSOS = True
                # This updates the referenced_variable_ids, just in case
                # there is a variable that only appears in an
                # SOSConstraint, in which case this needs to be known
                # before we write the "bounds" section (Cplex does not
                # handle this correctly, Gurobi does)
                self.printSOS(symbol_map,
                              labeler,
                              variable_symbol_map,
                              soscondata,
                              SOSlines)

        #
        # Bounds
        #

        output.append("bounds\n")

        # Scan all variables even if we're only writing a subset of them.
        # required because we don't store maps by variable type currently.

        # Track the number of integer and binary variables, so you can
        # output their status later.
        integer_vars = []
        binary_vars = []
        for vardata in variable_list:

            # TODO: We could just loop over the set of items in
            #       self._referenced_variable_ids, except this is
            #       a dictionary that is hashed by id(vardata)
            #       which would make the bounds section
            #       nondeterministic (bad for unit testing)
            if (not config.include_all_variable_bounds) and \
               (id(vardata) not in self._referenced_variable_ids):
                continue

            name_to_output = variable_symbol_dictionary[id(vardata)]
            if name_to_output == "e":
                raise ValueError(
                    "Attempting to write variable with name 'e' in a CPLEX LP "
                    "formatted file will cause a parse failure due to confusion with "
                    "numeric values expressed in scientific notation")

            # track the number of integer and binary variables, so we know whether
            # to output the general / binary sections below.
            if vardata.is_binary():
                binary_vars.append(name_to_output)
            elif vardata.is_integer():
                integer_vars.append(name_to_output)
            elif not vardata.is_continuous():
                raise TypeError("Invalid domain type for variable with name '%s'. "
                                "Variable is not continuous, integer, or binary."
                                % (vardata.name))

            if vardata.fixed:
                if not config.output_fixed_variable_bounds:
                    raise ValueError(
                        "Encountered a fixed variable (%s) inside an active "
                        "objective or constraint expression on model %s, which is "
                        "usually indicative of a preprocessing error. Use the "
                        "IO-option 'output_fixed_variable_bounds=True' to suppress "
                        "this error and fix the variable by overwriting its bounds "
                        "in the LP file." % (vardata.name, model.name))
                if vardata.value is None:
                    raise ValueError("Variable cannot be fixed to a value of None.")
                vardata_lb = value(vardata.value)
                vardata_ub = value(vardata.value)

                output.append("   ")
                output.append(lb_string_template
                                      % (_no_negative_zero(vardata_lb)))
                output.append(name_to_output)
                output.append(ub_string_template
                                      % (_no_negative_zero(vardata_ub)))
            else:
                vardata_lb = _get_bound(vardata.lb)
                vardata_ub = _get_bound(vardata.ub)

                # Pyomo assumes that the default variable bounds are -inf and +inf
                output.append("   ")
                if vardata.has_lb():
                    output.append(lb_string_template
                                      % (_no_negative_zero(vardata_lb)))
                else:
                    output.append(" -inf <= ")

                output.append(name_to_output)
                if vardata.has_ub():
                    output.append(ub_string_template
                                      % (_no_negative_zero(vardata_ub)))
                else:
                    output.append(" <= +inf\n")

        if len(integer_vars) > 0:

            output.append("general\n")
            for var_name in integer_vars:
                output.append('  %s\n' % var_name)

        if len(binary_vars) > 0:

            output.append("binary\n")
            for var_name in binary_vars:
                output.append('  %s\n' % var_name)


        # Write the SOS section
        output.append( "".join(SOSlines) )

        #
        # wrap-up
        #
        output.append("end\n")
        output_file.write( "".join(output) )

        # Clean up the symbol map to only contain variables referenced
        # in the active constraints **Note**: warm start method may
        # rely on this for choosing the set of potential warm start
        # variables
        vars_to_delete = set(variable_symbol_map.byObject.keys()) - \
                         set(self._referenced_variable_ids.keys())
        sm_byObject = symbol_map.byObject
        sm_bySymbol = symbol_map.bySymbol
        var_sm_byObject = variable_symbol_map.byObject
        for varid in vars_to_delete:
            symbol = var_sm_byObject[varid]
            del sm_byObject[varid]
            del sm_bySymbol[symbol]
        del variable_symbol_map

        return symbol_map

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################

class _objectOrderer(object):
    def __init__(self, initial_order):
        self.objOrder = {id(_): i for i,_ in enumerate(initial_order)}
        self.n = len(initial_order) - 1

    def object_sort_key(self, term):
        _id = id(term)
        if _id in self.objOrder:
            return self.objOrder[_id]
        else:
            self.n += 1
            n1 = self.objOrder[_id] = self.n
            return n1

    def linear_sort_key(self, term):
        _id = id(term[1])
        if _id in self.objOrder:
            return self.objOrder[_id]
        else:
            self.n += 1
            n1 = self.objOrder[_id] = self.n
            return n1

    def quadratic_sort_key(self, term):
        return self.objOrder[term[1]], self.objOrder[term[2]],

    def generate_quad_term(self, c, term, is_objective):
        _id = id(term[0])
        if _id in self.objOrder:
            n0 = self.objOrder[_id]
        else:
            self.n += 1
            n0 = self.objOrder[_id] = self.n

        _id = id(term[1])
        if _id in self.objOrder:
            n1 = self.objOrder[_id]
        else:
            self.n += 1
            n1 = self.objOrder[_id] = self.n

        # See note in _print_standard_repn for discussion on the
        # coefficient doubling for objectives
        if n1 < n0:
            return c*2. if is_objective else c, term[1], term[0]
        else:
            return c*2. if is_objective else c, term[0], term[1]


@WriterFactory.register('new_lp2', 'Generate the corresponding CPLEX LP file')
class ProblemWriter_cpxlp(object):

    CONFIG = ConfigBlock()
    CONFIG.declare('symbolic_solver_labels', ConfigValue(
        default=False,
        domain=bool,
        doc="""Use full Pyomo component names in the LP file rather than
        shortened symbols (slower, but useful for debugging)."""
    ))
    CONFIG.declare('labeler', ConfigValue(
        default=None,
    ))
    CONFIG.declare('skip_trivial_constraints', ConfigValue(
        default=False,
        domain=bool,
        doc="""Skip writing constraints whose body section is fixed (i.e., no
        variables)"""
    ))
    CONFIG.declare('output_fixed_variable_bounds', ConfigValue(
        default=False,
        domain=bool,
    ))
    CONFIG.declare('include_all_variable_bounds', ConfigValue(
        default=False,
        domain=bool,
        doc="""If False, unused variables will not be included in the
        LP file. Otherwise, include all variables in the bounds sections."""
    ))
    CONFIG.declare('file_determinism', ConfigValue(
        default=1,
        domain=In([0,1,2]),
        doc="""Specify how much effort do we want to put into ensuring the
        LP file is written deterministically for a Pyomo model:
          - 0 : None
          - 1 : sort keys of indexed components (default)
          - 2 : sort keys AND sort names (over declaration order)
        """
    ))
    CONFIG.declare('row_order', ConfigValue(
        default=None,
        domain=list,
        doc="Specify the order that constraints appear in the output"
    ))
    CONFIG.declare('column_order', ConfigValue(
        default=None,
        domain=list,
        doc="Specify the order that variables appear in the output"
    ))
    CONFIG.declare('force_objective_constant', ConfigValue(
        default=False,
        domain=bool,
        doc="""Make sure the ONE_VAR_CONSTANT variable appears in the
        objective even if the constant part of the objective is zero"""
    ))
    CONFIG.declare('allow_quadratic_objective', ConfigValue(
        default=True,
        domain=bool,
    ))
    CONFIG.declare('allow_quadratic_constraints', ConfigValue(
        default=True,
        domain=bool,
    ))
    CONFIG.declare('allow_sos1', ConfigValue(
        default=True,
        domain=bool,
    ))
    CONFIG.declare('allow_sos2', ConfigValue(
        default=True,
        domain=bool,
    ))

    def __init__(self):

        super(ProblemWriter_cpxlp, self).__init__()


        # Per ticket #4319, we are using %.17g, which mocks the
        # behavior of using %r (i.e., float('%r'%<number>) == <number>) with
        # the added benefit of outputting (+/-). The only case where this
        # fails to mock the behavior of %r is for large (long) integers (L),
        # which is a rare case to run into and is probably indicative of
        # other issues with the model.
        # *** NOTE ***: If you use 'r' or 's' here, it will break code that
        #               relies on using '%+' before the formatting character
        #               and you will need to go add extra logic to output
        #               the number's sign.
        self._precision_string = '.17g'
        self.linear_term_template =    '%+'+self._precision_string+' %s\n'
        self.squared_term_template =   '%+'+self._precision_string+' %s ^ 2\n'
        self.quadratic_term_template = '%+'+self._precision_string+' %s * %s\n'

        self.eq_string_template = "= %"+self._precision_string+'\n'

        self.obj_string_template =         '%+'+self._precision_string+' %s\n'
        self.sos_template_string =         "%s:%"+self._precision_string+"\n"
        self.geq_string_template =         ">= %"+self._precision_string+'\n\n'
        self.leq_string_template =         "<= %"+self._precision_string+'\n\n'
        self.lb_string_template =          "%"+self._precision_string+" <= "
        self.ub_string_template =          " <= %"+self._precision_string+"\n"

    def __call__(self,
                 model,
                 output_filename,
                 io_options=None,
                 **kwds):

        config = self.CONFIG(io_options)(kwds)

        if config.symbolic_solver_labels and (config.labeler is not None):
            raise ValueError("ProblemWriter_cpxlp: Using both the "
                             "'symbolic_solver_labels' and 'labeler' "
                             "I/O options is forbidden")

        #
        # Create labeler
        #
        if config.symbolic_solver_labels:
            labeler = TextLabeler()
        elif config.labeler is not None:
            labeler = config.labeler
        else:
            labeler = NumericLabeler('x')

        if output_filename is None:
            output_filename = model.name + ".lp"

        # when sorting, there are a non-trivial number of temporary
        # objects created. these all yield non-circular references, so
        # disable GC - the overhead is non-trivial, and because
        # references are non-circular, in CPython everything will be
        # collected immediately anyway.
        with PauseGC() as pgc:
            with open(output_filename, "w") as output:
                symbol_map = self._print_model_LP(
                    config, model, output, labeler)

        return output_filename, symbol_map

    def _print_model_LP(self, config, model, output, labeler):
        eq_template = self.eq_string_template
        geq_template = self.geq_string_template
        leq_template = self.leq_string_template

        symbol_map = SymbolMap()

        # Prep the variable ordering dictionary
        varOrder = _objectOrderer(config.column_order)
        conOrder = _objectOrderer(config.row_order)

        objective_generator = model.component_data_objects(
            Objective, active=True)
        obj = next(objective_generator, None)
        if obj is None:
            raise ValueError(
                "ERROR: No objectives defined for input model. "
                "Cannot write legal LP file.")
        if next(objective_generator, None) is not None:
            raise ValueError(
                "More than one active objective defined for input "
                "model '%s'; Cannot write legal LP file\n"
                "Objectives: %s" % (model.name, ' '.join(
                    o.name for o in model.component_data_objects(
                        Objective, active=True))))

        repn = generate_standard_repn(obj.expr)
        if repn.nonlinear_expr is not None:
            raise RuntimeError(
                "Cannot write legal LP file.  Objective '%s' "
                "has nonlinear terms that are not quadratic."
                % obj.name)
        if any(repn.quadratic_coefs):
            if not config.allow_quadratic_objective:
                raise RuntimeError(
                    "Objective contains quadratic terms and "
                    "allow_quadratic_objective is False.  Cannot write "
                    "the requested LP file.\n\tObjective at issue: %s."
                    % (obj.name,))
        elif not any(repn.linear_coefs):
            logger.warning("Constant objective detected, replacing "
                           "with a placeholder to prevent solver failure.")
            config.force_objective_constant = True

        if obj.is_minimizing():
            output.write("min \n%s:\n" % symbol_map.getSymbol(obj, labeler))
        else:
            output.write("max \n%s:\n" % symbol_map.getSymbol(obj, labeler))
        self._print_standard_repn(
            output, repn, True, varOrder, symbol_map, labeler)

        # Constraints
        #
        # If there are no non-trivial constraints, you'll end up with an empty
        # constraint block. CPLEX is OK with this, but GLPK isn't. And
        # eliminating the constraint block (i.e., the "s.t." line) causes GLPK
        # to whine elsewhere. Output a warning if the constraint block is empty,
        # so users can quickly determine the cause of the solve failure.

        output.write("\ns.t.\n\n")
        for con in sorted(model.component_data_objects(Constraint, active=True),
                          key=conOrder.object_sort_key):
            repn = generate_standard_repn(con.body)
            if repn.nonlinear_expr is not None:
                raise RuntimeError(
                    "Cannot write legal LP file.  Constraint '%s' "
                    "has nonlinear terms that are not quadratic."
                    % con.name)
            if any(repn.quadratic_coefs):
                if not config.allow_quadratic_constraints:
                    raise RuntimeError(
                        "Constraint contains quadratic terms and "
                        "allow_quadratic_constraints is False.  Cannot write "
                        "the requested LP file.\n\nConstraint at issue: %s."
                    % (con.name,))
            elif not any(repn.linear_coefs):
                if config.skip_trivial_constraints:
                    continue

            if con.equality:
                assert value(con.lower) == value(con.upper)
                output.write(sym(con, labeler, 'c_e_%s_') + ":\n")
                self._print_standard_repn(
                    output, repn, False, varOrder, symbol_map, labeler)
                output.write(eq_template
                             % (_get_bound(con.lower) - repn.constant))
            else:
                lb = _get_bound(con.lower)
                ub = _get_bound(con.upper)
                if lb is not None:
                    if ub is None:
                        output.write(sym(con, labeler, 'c_l_%s_') + ":\n")
                    else:
                        output.write(sym(con, labeler, 'r_l_%s_') + ":\n")
                    self._print_standard_repn(
                        output, repn, False, varOrder, symbol_map, labeler)
                    output.write(geq_template % (lb - repn.constant))
                if ub is not None:
                    if lb is None:
                        output.write(sym(con, labeler, 'c_u_%s_') + ":\n")
                    else:
                        output.write(sym(con, labeler, 'r_u_%s_') + ":\n")
                    self._print_standard_repn(
                        output, repn, False, varOrder, symbol_map, labeler)
                    output.write(leq_template % (ub - repn.constant))

        # the CPLEX LP format doesn't allow constants in the objective (or
        # constraint body), which is a bit silly.  To avoid painful
        # book-keeping, we introduce the following "variable", constrained
        # to the value 1.  This is used when quadratic terms are present.
        # worst-case, if not used, is that CPLEX easily pre-processes it out.
        output.write('c_e_ONE_VAR_CONSTANT_: \n'
                     'ONE_VAR_CONSTANT = 1.0\n\n')


    def _print_standard_repn(self, output, repn, is_objective,
                             varOrder, symbol_map, labeler):
        sym = symbol_map.getSymbol
        lvars = repn.linear_vars
        lcoef = repn.linear_coefs
        qvars = repn.quadratic_vars
        qcoef = repn.quadratic_coefs

        linear_term_template = self.linear_term_template
        squared_term_template = self.squared_term_template
        quadratic_term_template = self.quadratic_term_template

        constant_expression = True

        # Write the linear coeffs:
        if any(lcoef):
            constant_expression = False
            output.write(''.join(
                linear_term_template % (term[0], sym(term[1], labeler))
                for term in sorted(
                    ( (c, lvars[i]) for i,c in enumerate(lcoef) if c != 0 ),
                    key=varOrder.linear_sort_key )
            ))

        # Write the quadratic coefficients
        #
        # Note that the LP format appears to require that all quadratic
        # terms are doubled and then the entire quadratic section is
        # divided by 2, but ONLY for the objective expression.  This was
        # discovered through trial-and-error, plus the help of ILog
        # CPlex 8.0 User's Manual, p197.  The coefficient doubling is
        # handled by generate_quad_term().
        #
        if any(qcoef):
            constant_expression = False
            output.write("+ [\n" + ''.join(
                # Special handling if the two variables in the quadratic
                # term are the same
                squared_term_template % (term[0], sym(term[1], labeler))
                if term[1] is term[2] else
                quadratic_term_template % (
                    term[0], sym(term[1], labeler), sym(term[2], labeler))

                for term in sorted(
                        varOrder.generate_quad_term(c, qvars[i], is_objective)
                        for i,c in enumerate(qcoef) if c != 0,
                    key=varOrder.quadratic_sort_key )
                ) + ("] / 2\n" if is_objective else "]\n"))

        if is_objective:
            # Currently, it appears that we only need to print the
            # constant offset term multiplied by a variable for
            # objectives.
            if repn.constant or config.force_objective_constant:
                output.write(linear_term_template % (
                    repn.constant, 'ONE_VAR_CONSTANT' ))
        elif constant_expression:
            # If we made it to here we are outputing trivial constraints
            # place 0 * ONE_VAR_CONSTANT on this side of the constraint
            # for the benefit of solvers like Glpk that cannot parse an
            # LP file without a variable on the left hand side.
            output.append(linear_term_template % (0, 'ONE_VAR_CONSTANT'))



    def printSOS(self,
                 symbol_map,
                 labeler,
                 variable_symbol_map,
                 soscondata,
                 output):
        """
        Prints the SOS constraint associated with the _SOSConstraintData object
        """
        sos_template_string = self.sos_template_string

        if hasattr(soscondata, 'get_items'):
            sos_items = list(soscondata.get_items())
        else:
            sos_items = list(soscondata.items())

        if len(sos_items) == 0:
            return

        level = soscondata.level

        output.append('%s: S%s::\n'
                          % (symbol_map.getSymbol(soscondata,labeler), level))

        for vardata, weight in sos_items:
            weight = _get_bound(weight)
            if weight < 0:
                raise ValueError(
                    "Cannot use negative weight %f "
                    "for variable %s is special ordered "
                    "set %s " % (weight, vardata.name, soscondata.name))
            if vardata.fixed:
                raise RuntimeError(
                    "SOSConstraint '%s' includes a fixed variable '%s'. This is "
                    "currently not supported. Deactive this constraint in order to "
                    "proceed." % (soscondata.name, vardata.name))
            self._referenced_variable_ids[id(vardata)] = vardata
            output.append(sos_template_string
                              % (variable_symbol_map.getSymbol(vardata),
                                 weight))

    def XXX_print_model_LP(self,
                        config,
                        model,
                        output_file,
                        labeler):

        eq_string_template = self.eq_string_template
        leq_string_template = self.leq_string_template
        geq_string_template = self.geq_string_template
        ub_string_template = self.ub_string_template
        lb_string_template = self.lb_string_template

        symbol_map = SymbolMap()
        variable_symbol_map = SymbolMap()
        # NOTE: we use createSymbol instead of getSymbol because we
        #       know whether or not the symbol exists, and don't want
        #       to the overhead of error/duplicate checking.
        # cache frequently called functions
        create_symbol_func = SymbolMap.createSymbol
        create_symbols_func = SymbolMap.createSymbols
        alias_symbol_func = SymbolMap.alias
        variable_label_pairs = []

        # populate the symbol map in a single pass.
        #objective_list, constraint_list, sosconstraint_list, variable_list \
        #    = self._populate_symbol_map(model,
        #                                symbol_map,
        #                                labeler,
        #                                variable_symbol_map,
        #                                file_determinism=file_determinism)
        sortOrder = SortComponents.unsorted
        if config.file_determinism >= 1:
            sortOrder = sortOrder | SortComponents.indices
            if config.file_determinism >= 2:
                sortOrder = sortOrder | SortComponents.alphabetical

        #
        # Create variable symbols (and cache the block list)
        #
        all_blocks = []
        variable_list = []
        #
        # WEH - TODO:  See if this is faster
        # NOTE: This loop doesn't find all of the variables.  :(
        #
        #for block in model.block_data_objects(active=True,
        #                                      sort=sortOrder):
        #
        #    all_blocks.append(block)
        #
        #    for vardata in block.component_data_objects(
        #            Var,
        #            active=True,
        #            sort=sortOrder,
        #            descend_into=False):
        #
        #        variable_list.append(vardata)
        #        variable_label_pairs.append(
        #            (vardata,create_symbol_func(symbol_map,
        #                                        vardata,
        #                                        labeler)))
        all_blocks = list( model.block_data_objects(
                active=True, sort=sortOrder) )
        variable_list = list( model.component_data_objects(
                Var, sort=sortOrder) )
        variable_label_pairs = list(
            (vardata, create_symbol_func(symbol_map, vardata, labeler))
            for vardata in variable_list )
        variable_symbol_map.addSymbols(variable_label_pairs)

        # and extract the information we'll need for rapid labeling.
        object_symbol_dictionary = symbol_map.byObject
        variable_symbol_dictionary = variable_symbol_map.byObject

        # cache - these are called all the time.
        print_expr_canonical = self._print_expr_canonical

        # print the model name and the source, so we know roughly where
        # it came from.
        #
        # NOTE: this *must* use the "\* ... *\" comment format: the GLPK
        # LP parser does not correctly handle other formats (notably, "%").
        output = []
        output.append(
            "\\* Source Pyomo model name=%s *\\\n\n" % (model.name,) )

        #
        # Objective
        #

        supports_quadratic_objective = config.allow_quadratic_objective

        numObj = 0
        onames = []
        for block in all_blocks:

            gen_obj_repn = getattr(block, "_gen_obj_repn", True)

            # Get/Create the ComponentMap for the repn
            if not hasattr(block,'_repn'):
                block._repn = ComponentMap()
            block_repn = block._repn

            for objective_data in block.component_data_objects(
                    Objective,
                    active=True,
                    sort=sortOrder,
                    descend_into=False):

                numObj += 1
                onames.append(objective_data.name)
                if numObj > 1:
                    raise ValueError(
                        "More than one active objective defined for input "
                        "model '%s'; Cannot write legal LP file\n"
                        "Objectives: %s" % (model.name, ' '.join(onames)))

                create_symbol_func(symbol_map,
                                   objective_data,
                                   labeler)

                symbol_map.alias(objective_data, '__default_objective__')
                if objective_data.is_minimizing():
                    output.append("min \n")
                else:
                    output.append("max \n")

                if gen_obj_repn:
                    repn = generate_standard_repn(objective_data.expr)
                    block_repn[objective_data] = repn
                else:
                    repn = block_repn[objective_data]

                degree = repn.polynomial_degree()

                if degree == 0:
                    logger.warning("Constant objective detected, replacing "
                          "with a placeholder to prevent solver failure.")
                    config.force_objective_constant = True
                elif degree == 2:
                    if not supports_quadratic_objective:
                        raise RuntimeError(
                            "Selected solver is unable to handle "
                            "objective functions with quadratic terms. "
                            "Objective at issue: %s."
                            % objective_data.name)
                elif degree is None:
                    raise RuntimeError(
                        "Cannot write legal LP file.  Objective '%s' "
                        "has nonlinear terms that are not quadratic."
                        % objective_data.name)

                output.append(
                    object_symbol_dictionary[id(objective_data)]+':\n')

                offset = print_expr_canonical(
                    repn,
                    output,
                    object_symbol_dictionary,
                    variable_symbol_dictionary,
                    True,
                    config)

        if numObj == 0:
            raise ValueError(
                "ERROR: No objectives defined for input model. "
                "Cannot write legal LP file.")

        # Constraints
        #
        # If there are no non-trivial constraints, you'll end up with an empty
        # constraint block. CPLEX is OK with this, but GLPK isn't. And
        # eliminating the constraint block (i.e., the "s.t." line) causes GLPK
        # to whine elsewhere. Output a warning if the constraint block is empty,
        # so users can quickly determine the cause of the solve failure.

        output.append("\n")
        output.append("s.t.\n")
        output.append("\n")

        have_nontrivial = False

        supports_quadratic_constraint = config.allow_quadratic_constraints

        def constraint_generator():
            for block in all_blocks:

                gen_con_repn = getattr(block, "_gen_con_repn", True)

                # Get/Create the ComponentMap for the repn
                if not hasattr(block,'_repn'):
                    block._repn = ComponentMap()
                block_repn = block._repn

                for constraint_data in block.component_data_objects(
                        Constraint,
                        active=True,
                        sort=sortOrder,
                        descend_into=False):

                    if (not constraint_data.has_lb()) and \
                       (not constraint_data.has_ub()):
                        assert not constraint_data.equality
                        continue # non-binding, so skip

                    if constraint_data._linear_canonical_form:
                        repn = constraint_data.canonical_form()
                    elif gen_con_repn:
                        repn = generate_standard_repn(constraint_data.body)
                        block_repn[constraint_data] = repn
                    else:
                        repn = block_repn[constraint_data]

                    yield constraint_data, repn

        if config.row_order is not None:
            sorted_constraint_list = list(constraint_generator())
            sorted_constraint_list.sort(key=lambda x: config.row_order[x[0]])
            def yield_all_constraints():
                for data, repn in sorted_constraint_list:
                    yield data, repn
        else:
            yield_all_constraints = constraint_generator

        # FIXME: This is a hack to get nested blocks working...
        for constraint_data, repn in yield_all_constraints():
            have_nontrivial = True

            degree = repn.polynomial_degree()

            #
            # Write constraint
            #

            # There are conditions, e.g., when fixing variables, under which
            # a constraint block might be empty.  Ignore these, for both
            # practical reasons and the fact that the CPLEX LP format
            # requires a variable in the constraint body.  It is also
            # possible that the body of the constraint consists of only a
            # constant, in which case the "variable" of
            if degree == 0:
                if config.skip_trivial_constraints:
                    continue
            elif degree == 2:
                if not supports_quadratic_constraint:
                    raise ValueError(
                        "Solver unable to handle quadratic expressions. Constraint"
                        " at issue: '%s'" % (constraint_data.name))
            elif degree is None:
                raise ValueError(
                    "Cannot write legal LP file.  Constraint '%s' has a body "
                    "with nonlinear terms." % (constraint_data.name))

            # Create symbol
            con_symbol = create_symbol_func(symbol_map, constraint_data, labeler)

            if constraint_data.equality:
                assert value(constraint_data.lower) == \
                    value(constraint_data.upper)
                label = 'c_e_%s_' % con_symbol
                alias_symbol_func(symbol_map, constraint_data, label)
                output.append(label)
                output.append(':\n')
                offset = print_expr_canonical(repn,
                                              output,
                                              object_symbol_dictionary,
                                              variable_symbol_dictionary,
                                              False,
                                              config)
                bound = constraint_data.lower
                bound = _get_bound(bound) - offset
                output.append(eq_string_template
                                  % (_no_negative_zero(bound)))
                output.append("\n")
            else:
                if constraint_data.has_lb():
                    if constraint_data.has_ub():
                        label = 'r_l_%s_' % con_symbol
                    else:
                        label = 'c_l_%s_' % con_symbol
                    alias_symbol_func(symbol_map, constraint_data, label)
                    output.append(label)
                    output.append(':\n')
                    offset = print_expr_canonical(repn,
                                                  output,
                                                  object_symbol_dictionary,
                                                  variable_symbol_dictionary,
                                                  False,
                                                  config)
                    bound = constraint_data.lower
                    bound = _get_bound(bound) - offset
                    output.append(geq_string_template
                                      % (_no_negative_zero(bound)))
                else:
                    assert constraint_data.has_ub()

                if constraint_data.has_ub():
                    if constraint_data.has_lb():
                        label = 'r_u_%s_' % con_symbol
                    else:
                        label = 'c_u_%s_' % con_symbol
                    alias_symbol_func(symbol_map, constraint_data, label)
                    output.append(label)
                    output.append(':\n')
                    offset = print_expr_canonical(repn,
                                                  output,
                                                  object_symbol_dictionary,
                                                  variable_symbol_dictionary,
                                                  False,
                                                  config)
                    bound = constraint_data.upper
                    bound = _get_bound(bound) - offset
                    output.append(leq_string_template
                                      % (_no_negative_zero(bound)))
                else:
                    assert constraint_data.has_lb()

            # A simple hack to avoid caching super large files
            if len(output) > 1024:
                output_file.write( "".join(output) )
                output = []

        if not have_nontrivial:
            logger.warning('Empty constraint block written in LP format '  \
                  '- solver may error')

        # the CPLEX LP format doesn't allow constants in the objective (or
        # constraint body), which is a bit silly.  To avoid painful
        # book-keeping, we introduce the following "variable", constrained
        # to the value 1.  This is used when quadratic terms are present.
        # worst-case, if not used, is that CPLEX easily pre-processes it out.
        prefix = ""
        output.append('%sc_e_ONE_VAR_CONSTANT: \n' % prefix)
        output.append('%sONE_VAR_CONSTANT = 1.0\n' % prefix)
        output.append("\n")

        # SOS constraints
        #
        # For now, we write out SOS1 and SOS2 constraints in the cplex format
        #
        # All Component objects are stored in model._component, which is a
        # dictionary of {class: {objName: object}}.
        #
        # Consider the variable X,
        #
        #   model.X = Var(...)
        #
        # We print X to CPLEX format as X(i,j,k,...) where i, j, k, ... are the
        # indices of X.
        #
        SOSlines = []
        sos1 = config.allow_sos1
        sos2 = config.allow_sos2
        writtenSOS = False
        for block in all_blocks:

            for soscondata in block.component_data_objects(
                    SOSConstraint,
                    active=True,
                    sort=sortOrder,
                    descend_into=False):

                create_symbol_func(symbol_map, soscondata, labeler)

                level = soscondata.level
                if (level == 1 and not sos1) or \
                   (level == 2 and not sos2) or \
                   (level > 2):
                    raise ValueError(
                        "Solver does not support SOS level %s constraints" % (level))
                if writtenSOS == False:
                    SOSlines.append("SOS\n")
                    writtenSOS = True
                # This updates the referenced_variable_ids, just in case
                # there is a variable that only appears in an
                # SOSConstraint, in which case this needs to be known
                # before we write the "bounds" section (Cplex does not
                # handle this correctly, Gurobi does)
                self.printSOS(symbol_map,
                              labeler,
                              variable_symbol_map,
                              soscondata,
                              SOSlines)

        #
        # Bounds
        #

        output.append("bounds\n")

        # Scan all variables even if we're only writing a subset of them.
        # required because we don't store maps by variable type currently.

        # Track the number of integer and binary variables, so you can
        # output their status later.
        integer_vars = []
        binary_vars = []
        for vardata in variable_list:

            # TODO: We could just loop over the set of items in
            #       self._referenced_variable_ids, except this is
            #       a dictionary that is hashed by id(vardata)
            #       which would make the bounds section
            #       nondeterministic (bad for unit testing)
            if (not config.include_all_variable_bounds) and \
               (id(vardata) not in self._referenced_variable_ids):
                continue

            name_to_output = variable_symbol_dictionary[id(vardata)]
            if name_to_output == "e":
                raise ValueError(
                    "Attempting to write variable with name 'e' in a CPLEX LP "
                    "formatted file will cause a parse failure due to confusion with "
                    "numeric values expressed in scientific notation")

            # track the number of integer and binary variables, so we know whether
            # to output the general / binary sections below.
            if vardata.is_binary():
                binary_vars.append(name_to_output)
            elif vardata.is_integer():
                integer_vars.append(name_to_output)
            elif not vardata.is_continuous():
                raise TypeError("Invalid domain type for variable with name '%s'. "
                                "Variable is not continuous, integer, or binary."
                                % (vardata.name))

            if vardata.fixed:
                if not config.output_fixed_variable_bounds:
                    raise ValueError(
                        "Encountered a fixed variable (%s) inside an active "
                        "objective or constraint expression on model %s, which is "
                        "usually indicative of a preprocessing error. Use the "
                        "IO-option 'output_fixed_variable_bounds=True' to suppress "
                        "this error and fix the variable by overwriting its bounds "
                        "in the LP file." % (vardata.name, model.name))
                if vardata.value is None:
                    raise ValueError("Variable cannot be fixed to a value of None.")
                vardata_lb = value(vardata.value)
                vardata_ub = value(vardata.value)

                output.append("   ")
                output.append(lb_string_template
                                      % (_no_negative_zero(vardata_lb)))
                output.append(name_to_output)
                output.append(ub_string_template
                                      % (_no_negative_zero(vardata_ub)))
            else:
                vardata_lb = _get_bound(vardata.lb)
                vardata_ub = _get_bound(vardata.ub)

                # Pyomo assumes that the default variable bounds are -inf and +inf
                output.append("   ")
                if vardata.has_lb():
                    output.append(lb_string_template
                                      % (_no_negative_zero(vardata_lb)))
                else:
                    output.append(" -inf <= ")

                output.append(name_to_output)
                if vardata.has_ub():
                    output.append(ub_string_template
                                      % (_no_negative_zero(vardata_ub)))
                else:
                    output.append(" <= +inf\n")

        if len(integer_vars) > 0:

            output.append("general\n")
            for var_name in integer_vars:
                output.append('  %s\n' % var_name)

        if len(binary_vars) > 0:

            output.append("binary\n")
            for var_name in binary_vars:
                output.append('  %s\n' % var_name)


        # Write the SOS section
        output.append( "".join(SOSlines) )

        #
        # wrap-up
        #
        output.append("end\n")
        output_file.write( "".join(output) )

        # Clean up the symbol map to only contain variables referenced
        # in the active constraints **Note**: warm start method may
        # rely on this for choosing the set of potential warm start
        # variables
        vars_to_delete = set(variable_symbol_map.byObject.keys()) - \
                         set(self._referenced_variable_ids.keys())
        sm_byObject = symbol_map.byObject
        sm_bySymbol = symbol_map.bySymbol
        var_sm_byObject = variable_symbol_map.byObject
        for varid in vars_to_delete:
            symbol = var_sm_byObject[varid]
            del sm_byObject[varid]
            del sm_bySymbol[symbol]
        del variable_symbol_map

        return symbol_map

