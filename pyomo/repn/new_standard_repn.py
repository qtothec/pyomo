from pyomo.repn.standard_repn import StandardRepn

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

_CONSTANT = 0
_MONOMIAL = 2
_LINEAR = 4
_GENERAL = 8

def profile(x):
    return x

class GeneralStandardExpressionVisitor(StreamBasedExpressionVisitor_allCallbacks):
    linearExprPool = None

    @profile
    def initializeWalker(self, expr):
        if type(expr) in native_types:
            return False, (_CONSTANT, expr)
        elif not expr.is_expression_type():
            if expr.is_fixed():
                return False, (_CONSTANT, value(expr))
            else:
                return False, (_MONOMIAL, 1, expr)
        return True, expr

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

        if child_type is not LinearExpression:
            return True, None

        print "COPY Linear Expr!"
        # Because we are going to modify the LinearExpression in this
        # walker, we need to make a copy of the LinearExpression from
        # the original expression tree.
        linear = self._get_linear()
        idMap = {}
        zeros = set()
        next_pos = 0

        linear.constant = child.constant
        for c, v in zip(child.linear_coefs, child.linear_vars):
            c = value(c)
            if not c:
                continue
            elif v.is_fixed():
                linear.constant += c*value(v)
            else:
                _id = id(v)
                if _id in idMap:
                    i = idMap[_id]
                    linear.linear_coefs[i] += c
                    if not linear.linear_coefs[i]:
                        zeros.add(i)
                else:
                    idMap[_id] = next_pos
                    next_pos += 1
                    linear.linear_vars.append(v)
                    linear.linear_coefs.append(c)
        if zeros:
            self._finalize_linear(zeros, linear)
        return False, (_LINEAR, linear)

    def afterChild(self, node, child):
        pass

    @profile
    def enterNode(self, node):
        if node._precedence() != SumExpressionBase.PRECEDENCE:
            return node._args_, []
        return node._args_, ({}, set(), self._get_linear())

    @profile
    def acceptChildResult(self, node, data, child_result):
        if data.__class__ is list:
            # General expression... cache the child result until the end
            data.append(child_result)
        else:
            # Linear Expression
            idMap, zeros, linear = data
            child_type = child_result[0]
            if child_type is _MONOMIAL:
                _id = id(child_result[2])
                if _id not in idMap:
                    idMap[_id] = idx = len(linear.linear_coefs)
                    linear.linear_coefs.append(child_result[1])
                    linear.linear_vars.append(child_result[2])
                else:
                    idx = idMap[_id]
                    linear.linear_coefs[idx] += child_result[1]
                if not linear.linear_coefs[idx]:
                    zeros.add(idx)
            elif child_type is _CONSTANT:
                linear.constant += child_result[1]
            elif child_type is _LINEAR:
                child = child_result[1]
                next_id = len(linear.linear_coefs)
                linear.constant += child.constant
                for i in xrange(len(child.linear_coefs)):
                    _id = id(child.linear_vars[i])
                    if _id not in idMap:
                        idMap[_id] = idx = next_id
                        next_id += 1
                        linear.linear_coefs.append(child.linear_coefs[i])
                        linear.linear_vars.append(child.linear_vars[i])
                    else:
                        idx = idMap[_id]
                        linear.linear_coefs[idx] += child.linear_coefs[i]
                    if not linear.linear_coefs[idx]:
                        zeros.add(idx)
            elif child_type is _GENERAL:
                print "ERROR"
                # General (nonlinear) sub-expression.  We will convert
                # this from a Linear expression to a general summation
                if zeros:
                    self._finalize_linear(zeros, linear)
                if linear.constant or linear.linear_vars:
                    if not linear.linear_coefs:
                        const, linear.constant = linear.constant, 0
                        self.linearExprPool = (self.linearExprPool, linear)
                        return [(_CONSTANT, const), child_result]
                    return [(_LINEAR, linear), child_result]
                else:
                    return [child_result]

        return data

    @profile
    def exitNode(self, node, data):
        if data.__class__ is tuple:
            # Linear expression
            _, zeros, linear = data
            if zeros:
                self._finalize_linear(zeros, linear)
                if not linear.linear_coefs:
                    const, linear.constant = linear.constant, 0
                    self.linearExprPool = (self.linearExprPool, linear)
                    return (_CONSTANT, const)
            return (_LINEAR, linear)
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
        print "exit general"
        args = tuple(_[1]*_[2] if _[0] is _MONOMIAL else _[1])
        if all(_[0] is _CONSTANT for _ in data):
            return node._apply_operation(args)
        return (_GENERAL, node.create_node_with_local_data(args))

    @profile
    def finalizeResult(self, result):
        result_type = result[0]
        ans = StandardRepn()
        if result_type is _LINEAR:
            expr = result[1]
            ans.constant = expr.constant
            ans.linear_vars = tuple(expr.linear_vars)
            ans.linear_coefs = tuple(expr.linear_coefs)
            expr.constant = 0
            expr.linear_vars = []
            expr.linear_coefs = []
            self.linearExprPool = (self.linearExprPool, expr)
        elif result_type is _GENERAL:
            print "TODO: Separate Linear and Nonlinear terms"
            if isinstance(expr, SumExpressionBase):
                linear_terms = []
                linear = self._get_linear()
                idMap = {}
                zeros = set()
                for i,term in enumerate(expr._args_):
                    term_type = type(term)
                    if term_type in native_types:
                        term_type = _CONSTANT
                    elif term_type is MonomialTermExpression:
                        term_type = _MONOMIAL
                    elif term_type is LinearExpression:
                        term_type = _LINEAR
                    else:
                        continue
                    if self._update_linear_expr(
                            idMap, zeros, linear, term_type, term):
                        linear_terms.append(i)
                if zeros:
                    self._finalize_linear(zeros, linear)
                ans.constant = linear.constant
                ans.linear_vars = tuple(linear.linear_vars)
                ans.linear_coefs = tuple(linear.linear_coefs)
                for i in reversed(linear_terms):
                    expr._args_.pop(i)
                expr._nargs = len(expr._args_)
            ans.nonlinear_expr = expr
            ans.nonlinear_vars = list(identify_variables(expr))
        elif result_type is _MONOMIAL:
            print "FINALIZE monomial"
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
            return LinearExpression()

    @profile
    def _finalize_linear(zeros, linear):
        print "NOTE! finalize linear"
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

class QuadraticStandardExpressionVisitor(GeneralStandardExpressionVisitor):
    pass
