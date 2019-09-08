from __future__ import division

from pyomo.common.errors import DeveloperError
from pyomo.core.expr import (
    native_types, value,
)
from pyomo.core.expr.current import (
    ProductExpression, MonomialTermExpression, SumExpressionBase,
    LinearExpression, identify_variables,
)
from pyomo.core.expr.visitor import (
    StreamBasedExpressionVisitor, StreamBasedExpressionVisitor_allCallbacks
)

try:
    from pyomo.core.expr.current import DivisionExpression
except ImportError:
    class DivisionExpression(object): pass

from pyomo.repn.standard_repn import StandardRepn


_CONSTANT = 0
_MONOMIAL = 2
_LINEAR = 4
_GENERAL = 8

if not hasattr(__builtins__, 'profile'):
    def profile(x):
        return x

REMOVE_ZERO_COEF = False
INLINE = True

class _linearRepn(object):
    __slots__ = ('coef','vars','const')

    def __init__(self):
        self.coef = {}
        self.vars = []
        self.const = 0

    def merge(self, other):
        self.const += other.const
        for v in other.vars:
            _id = id(v)
            coef = other.coef[_id]
            if _id in self.coef:
                self.coef[_id] += coef
            elif coef:
                self.coef[_id] = coef
                self.vars.append(v)
        other.const = 0
        other.coef.clear()
        other.vars = []

    def toLinearExpr(self):
        ans = LinearExpression()
        ans.constant, self.const = self.const, 0
        ans.linear_coefs = list(self.coef[id(v)] for v in self.vars)
        ans.linear_vars, self.vars = self.vars, []
        if REMOVE_ZERO_COEF:
            try:
                i = 0
                while 1:
                    i = ans.linear_coefs.index(0,i)
                    ans.linear_coefs.pop(i)
                    ans.linear_vars.pop(i)
            except ValueError:
                pass
        self.coef.clear()
        return ans

    def fromLinearExpr(self, linear):
        self.const += linear.constant
        for i,coef in enumerate(linear.linear_coefs):
            if not coef:
                continue
            v = linear.linear_vars[i]
            _id = id(v)
            if _id in linear.coef:
                linear.coef[_id] += coef
            else:
                linear.coef[_id] = coef
                linear.vars.append(v)


class GeneralStandardExpressionVisitor_streambased(
        StreamBasedExpressionVisitor_allCallbacks):
    linearExprPool = None

    @profile
    def initializeWalker(self, expr):
        walk, result = self.beforeChild(None, expr)
        if not walk:
            return False, self.finalizeResult(result)
        return True, None

    @profile
    def beforeChild(self, node, child):
        child_type = child.__class__
        if child_type in native_types:
            return False, (_CONSTANT, child)
        if not child.is_expression_type():
            if child.is_fixed():
                return False, (_CONSTANT, value(child))
            else:
                return False, (_MONOMIAL, 1, child)

        #
        # The following are performance optimizations for common
        # situations (Monomial terms and Product expressions that are in
        # fact monomial terms)
        #

        if child_type is MonomialTermExpression:
            arg1, arg2 = child._args_
            if arg1.__class__ not in native_types:
                arg1 = value(arg1)
            if arg2.is_fixed():
                return False, (_CONSTANT, child())
            else:
                return False, (_MONOMIAL, arg1, arg2)

        if child_type is LinearExpression:
            print("COPY Linear Expr!")
            # Because we are going to modify the LinearExpression in this
            # walker, we need to make a copy of the LinearExpression from
            # the original expression tree.
            linear = self._get_linear()
            linear.fromLinearExpr(child)
            return False, (_LINEAR, linear)

        # if child_type is ProductExpression:
        #     arg1, arg2 = child._args_
        #     if arg1.__class__ in native_types or (
        #             not arg1.is_expression_type() and arg1.is_fixed()):
        #         if arg2.__class__ in native_types:
        #             return False, (_CONSTANT, child())
        #         elif not arg2.is_expression_type():
        #             if arg2.is_fixed():
        #                 return False, (_CONSTANT, child())
        #             else:
        #                 return False, (_MONOMIAL, value(arg1), arg2)
        #     elif arg2.__class__ in native_types or (
        #             not arg2.is_expression_type() and arg2.is_fixed()):
        #         if not arg1.is_expression_type():
        #             # We know arg1 is not fixed...
        #             return False, (_MONOMIAL, value(arg2), arg1)

        return True, None

    def afterChild(self, node, child):
        pass

    @profile
    def enterNode(self, node):
        if isinstance(node, SumExpressionBase):
            return node._args_, self._get_linear()
        else:
            return node._args_, []

    @profile
    def acceptChildResult(self, node, data, child_result):
        if data.__class__ is list:
            # General expression... cache the child result until the end
            data.append(child_result)
        else:
            # Linear Expression
            child_type = child_result[0]
            if child_type is _MONOMIAL:
                _id = id(child_result[2])
                if _id in data.coef:
                    data.coef[_id] += child_result[1]
                else:
                    data.coef[_id] = child_result[1]
                    data.vars.append(child_result[2])
            elif child_type is _CONSTANT:
                data.const += child_result[1]
            elif child_type is _LINEAR:
                child = child_result[1]
                data.merge(child)
                self.linearExprPool = (self.linearExprPool, child)
            elif child_type is _GENERAL:
                if data.const or data.vars:
                    if not data.vars:
                        const, data.const = data.const, 0
                        self.linearExprPool = (self.linearExprPool, data)
                        return [(_CONSTANT, const), child_result]
                    return [(_LINEAR, data), child_result]
                else:
                    return [child_result]
        return data

    @profile
    def exitNode(self, node, data):
        if data.__class__ is _linearRepn:
            return (_LINEAR, data)
        #
        # General (nonlinear) expression...
        #
        # If all the arguments are constants, then evaluate this node
        # if all(_[0] is _CONSTANT for _ in data):
        #     return (
        #         _CONSTANT,
        #         node._apply_operation(tuple(_[1] for _ in data))
        #     )
        # We need special handling for Product/Division expressions
        if len(data) == 2:
            if isinstance(node, ProductExpression):
                if data[1][0] is _CONSTANT:
                    arg2, arg1 = data
                else:
                    arg1, arg2 = data
                if arg1[0] is _CONSTANT:
                    if arg2[0] is _MONOMIAL:
                        return (_MONOMIAL, arg1[1]*arg2[1], arg2[2])
                    elif arg2[0] is _LINEAR:
                        mul = arg1[1]
                        arg2[1].constant *= mul
                        for i in xrange(len(arg2[1].linear_coefs)):
                            arg2[1].linear_coefs[i] *= mul
                        return arg2
                    elif arg2[0] is _CONSTANT:
                        return (_CONSTANT, arg1[1]*arg2[1])
            elif isinstance(node, DivisionExpression):
                arg1, arg2 = data
                if arg2[0] is _CONSTANT:
                    div = arg2[1]
                    if arg1[0] is _MONOMIAL:
                        return (_MONOMIAL, (arg1[1]/div, arg1[2])
                        )
                    elif arg1[0] is _LINEAR:
                        arg1[1].constant /= div
                        for i in xrange(len(arg1[1].linear_coefs)):
                            arg1[1].linear_coefs[i] /= div
                        return arg1
                    elif arg1[0] is _CONSTANT:
                        return (_CONSTANT, arg1[1]/arg2[1])

        # We need to convert data to valid expression objects
        print("exit general")
        args = tuple( _[1]*_[2] if _[0] is _MONOMIAL
                      else _[1].toLinearExpr() if _[0] is _LINEAR
                      else _[1] for _ in data)
        if all(_[0] is _CONSTANT for _ in data):
            return node._apply_operation(args)
        return (_GENERAL, node.create_node_with_local_data(args))

    @profile
    def finalizeResult(self, result):
        result_type = result[0]
        ans = StandardRepn()
        if result_type is _LINEAR:
            expr = result[1]
            ans.constant, expr.const = expr.const, 0
            ans.linear_coefs = list(expr.coef[id(v)] for v in expr.vars)
            if REMOVE_ZERO_COEF:
                try:
                    i = 0
                    while 1:
                        i = ans.linear_coefs.index(0,i)
                        ans.linear_coefs.pop(i)
                        ans.linear_vars.pop(i)
                except ValueError:
                    pass
            expr.coef.clear()
            ans.linear_vars, expr.vars = expr.vars, []
            self.linearExprPool = (self.linearExprPool, expr)
        elif result_type is _GENERAL:
            print("TODO: Separate Linear and Nonlinear terms")
            expr = result[1]
            ans.nonlinear_expr = expr
            ans.nonlinear_vars = list(identify_variables(expr))
        elif result_type is _MONOMIAL:
            print("FINALIZE monomial")
            if result[1]:
                ans.linear_coefs = (result[1],)
                ans.linear_vars = (result[2],)
        elif result_type is _CONSTANT:
            ans.constant = result[1]
        else:
            raise DeveloperError("unknown result type")

        return ans

    @profile
    def _get_linear(self):
        if self.linearExprPool:
            self.linearExprPool, ans = self.linearExprPool
            return ans
        else:
            return _linearRepn()

    @profile
    def _finalize_linear(zeros, linear):
        print("NOTE! finalize linear")
        for i in reversed(sorted(zeros)):
            if not linear.linear_coefs[i]:
                linear.linear_coefs.pop(i)
                linear.linear_vars.pop(i)

    @profile
    def _update_linear_expr(
            self, idMap, zeros, linear, child_type, child):
        if child_type is _MONOMIAL:
            coef, child = child._args_
            child_type = self.VARIABLE
        else:
            coef = 1

        if child_type is self.VARIABLE:
            _id = id(child)
            if _id in idMap:
                linear.linear_coefs[idMap[_id]] += coef
            else:
                idMap[_id] = len(linear.linear_coefs)
                linear.linear_coefs.append(coef)
                linear.linear_vars.append(child)
        elif child_type is _CONSTANT:
            linear.constant += child
        elif child_type is _LINEAR:
            linear.constant += child.constant
            for c, v in zip(child.linear_coefs, child.linear_vars):
                if not c:
                    continue
                _id = id(v)
                if _id in idMap:
                    i = idMap[_id]
                    linear.linear_coefs[i] += c
                    if not linear.linear_coefs[i]:
                        zeros.add(i)
                else:
                    idMap[_id] = len(linear.linear_coefs)
                    linear.linear_coefs.append(c)
                    linear.linear_vars.append(v)
            child.constant = 0
            child.linear_vars = []
            child.linear_coefs = []
            self.linearExprPool = (self.linearExprPool, child)
        else:
            # Nonlinear expression
            return False
        return True


class GeneralStandardExpressionVisitor_inlined(object):
    linearExprPool = None

    @profile
    def walk_expression(self, expr):
        ##
        # Initialize walker
        ##
        expr_type = expr.__class__
        if expr_type in native_types:
            ans = StandardRepn()
            ans.constant = expr
            return ans
        elif not expr.is_expression_type():
            ans = StandardRepn()
            if expr.is_fixed():
                ans.constant = value(expr)
            else:
                ans.linear_vars = (expr,)
                ans.linear_coefs = (1,)
            return ans
        elif expr_type is MonomialTermExpression:
            ans = StandardRepn()
            coef = value(expr._args_[0])
            if coef:
                v = expr._args_[1]
                if v.is_fixed():
                    ans.constant = value(v) * coef
                else:
                    ans.linear_vars = (v,)
                    ans.linear_coefs = (coef,)
            return ans

        ##
        # Enter Node
        ##
        args = expr._args_
        if isinstance(expr, SumExpressionBase):
            data = self._get_linear()
        else:
            data = []
        #
        nargs = expr.nargs()
        node = expr
        child_idx = 0
        ptr = (None, node, args, nargs, data, child_idx)

        while 1:
            if child_idx < nargs:
                child = args[child_idx]
                child_idx += 1

                ##
                # Before Child
                ##
                child_type = child.__class__
                #
                # The following are performance optimizations for common
                # situations (Monomial terms)
                #
                if child_type is MonomialTermExpression:
                    arg1, arg2 = child._args_
                    if arg1.__class__ not in native_types:
                        arg1 = value(arg1)
                    if arg2.is_fixed():
                        node_result = (_CONSTANT, child())
                    else:
                        node_result = (_MONOMIAL, arg1, arg2)
                elif child_type in native_types:
                    node_result = (_CONSTANT, child)
                elif not child.is_expression_type():
                    if child.is_fixed():
                        node_result = (_CONSTANT, value(child))
                    else:
                        node_result = (_MONOMIAL, 1, child)
                elif child_type is LinearExpression:
                    print("COPY Linear Expr!")
                    # Because we are going to modify the
                    # LinearExpression in this walker, we need to make a
                    # copy of the LinearExpression from the original
                    # expression tree.
                    linear = self._get_linear()
                    linear.fromLinearExpr(child)
                    node_result = (_LINEAR, linear)
                else:
                    ##
                    # Descend into child
                    ##
                    ptr = ptr[:4] + (data, child_idx,)
                    ##
                    # Enter Node
                    ##
                    args = child._args_
                    if isinstance(child, SumExpressionBase):
                        data = self._get_linear()
                    else:
                        data = []
                    #
                    nargs = child.nargs()
                    node = child
                    child_idx = 0
                    ptr = (ptr, node, args, nargs, data, child_idx)
                    continue
            else:
                ##
                # Exit Node
                ##
                if data.__class__ is _linearRepn:
                    node_result = (_LINEAR, data)
                elif len(data) == 2:
                    node_result = None
                    if isinstance(node, ProductExpression):
                        if data[1][0] is _CONSTANT:
                            arg2, arg1 = data
                        else:
                            arg1, arg2 = data
                        if arg1[0] is _CONSTANT:
                            if arg2[0] is _MONOMIAL:
                                node_result = (
                                    _MONOMIAL, arg1[1]*arg2[1], arg2[2])
                            elif arg2[0] is _LINEAR:
                                mul = arg1[1]
                                arg2[1].constant *= mul
                                for i in xrange(len(arg2[1].linear_coefs)):
                                    arg2[1].linear_coefs[i] *= mul
                                node_result = arg2
                            elif arg2[0] is _CONSTANT:
                                node_result = (_CONSTANT, arg1[1]*arg2[1])
                    elif isinstance(node, DivisionExpression):
                        arg1, arg2 = data
                        if arg2[0] is _CONSTANT:
                            div = arg2[1]
                            if arg1[0] is _MONOMIAL:
                                node_result = (
                                    _MONOMIAL, (arg1[1]/div, arg1[2]))
                            elif arg1[0] is _LINEAR:
                                arg1[1].constant /= div
                                for i in xrange(len(arg1[1].linear_coefs)):
                                    arg1[1].linear_coefs[i] /= div
                                node_result = arg1
                            elif arg1[0] is _CONSTANT:
                                node_result = (_CONSTANT, arg1[1]/arg2[1])
                else:
                    node_result = None

                if node_result is None:
                    # We need to convert data to valid expression objects
                    print("exit general")
                    args = tuple( _[1]*_[2] if _[0] is _MONOMIAL
                                  else _[1].toLinearExpr() if _[0] is _LINEAR
                                  else _[1] for _ in data)
                    if all(_[0] is _CONSTANT for _ in data):
                        node_result = (_CONSTANT, node._apply_operation(args))
                    else:
                        node_result = (
                            _GENERAL, node.create_node_with_local_data(args))
                ##
                # Pop the node
                ##
                ptr = ptr[0]
                # If we have returned to the beginning, return the final
                # answer
                if ptr is None:
                    ##
                    # Finalize Result
                    ##
                    result_type = node_result[0]
                    ans = StandardRepn()
                    if result_type is _LINEAR:
                        expr = node_result[1]
                        ans.constant, expr.const = expr.const, 0
                        ans.linear_coefs = list(expr.coef[id(v)] for v in expr.vars)
                        if REMOVE_ZERO_COEF:
                            try:
                                i = 0
                                while 1:
                                    i = ans.linear_coefs.index(0,i)
                                    ans.linear_coefs.pop(i)
                                    ans.linear_vars.pop(i)
                            except ValueError:
                                pass
                        expr.coef.clear()
                        ans.linear_vars, expr.vars = expr.vars, []
                        self.linearExprPool = (self.linearExprPool, expr)
                    elif result_type is _GENERAL:
                        print("TODO: Separate Linear and Nonlinear terms")
                        expr = node_result[1]
                        ans.nonlinear_expr = expr
                        ans.nonlinear_vars = list(identify_variables(expr))
                    elif result_type is _MONOMIAL:
                        print("FINALIZE monomial")
                        if node_result[1]:
                            ans.linear_coefs = (node_result[1],)
                            ans.linear_vars = (node_result[2],)
                    elif result_type is _CONSTANT:
                        ans.constant = node_result[1]
                    else:
                        raise DeveloperError("unknown result type")
                    return ans

                # Not done yet, update node to point to the new active
                # node
                #child = node
                _, node, args, nargs, data, child_idx = ptr

            ##
            # Accept Child Result
            ##
            if data.__class__ is list:
                # General expression... cache the child result until the end
                data.append(node_result)
            else:
                # Linear Expression
                child_type = node_result[0]
                if child_type is _MONOMIAL:
                    _id = id(node_result[2])
                    if _id in data.coef:
                        data.coef[_id] += node_result[1]
                    else:
                        data.coef[_id] = node_result[1]
                        data.vars.append(node_result[2])
                elif child_type is _CONSTANT:
                    data.const += node_result[1]
                elif child_type is _LINEAR:
                    child = node_result[1]
                    data.merge(child)
                    self.linearExprPool = (self.linearExprPool, child)
                elif child_type is _GENERAL:
                    if data.const or data.vars:
                        if not data.vars:
                            const, data.const = data.const, 0
                            self.linearExprPool = (self.linearExprPool, data)
                            return [(_CONSTANT, const), node_result]
                        data = [(_LINEAR, data), node_result]
                    else:
                        data = [node_result]
            ##
            # After Child
            ##



    @profile
    def _get_linear(self):
        #if self.linearExprPool:
        try:
            self.linearExprPool, ans = self.linearExprPool
            return ans
        #else:
        except:
            return _linearRepn()


GeneralStandardExpressionVisitor = GeneralStandardExpressionVisitor_inlined \
    if INLINE else GeneralStandardExpressionVisitor_streambased

class QuadraticStandardExpressionVisitor(GeneralStandardExpressionVisitor):
    pass
