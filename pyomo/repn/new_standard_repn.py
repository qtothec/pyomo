from pyomo.repn.standard_repn import StandardRepn

from pyomo.core.expr import (
    native_types, value,
)
from pyomo.core.expr.current import (
    ProductExpression, MonomialTermExpression, SumExpressionBase,
    LinearExpression, identify_variables,
)
from pyomo.core.expr.visitor import (
    StreamBasedExpressionVisitor,
)

try:
    from pyomo.core.expr.current import DivisionExpression
except ImportError:
    class DivisionExpression(object): pass
    

class GeneralStandardExpressionVisitor(StreamBasedExpressionVisitor):
    linearExprPool = None
    monomialExprPool = None
    CONSTANT = 0
    VARIABLE = 1
    MONOMIAL = 2
    LINEAR = 4
    GENERAL = 8

    def initializeWalker(self, expr):
        if type(expr) in native_types:
            return False, (self.CONSTANT, expr)
        elif not expr.is_expression_type():
            if expr.is_fixed():
                return False, (self.CONSTANT, value(expr))
            else:
                return False, (self.VARIABLE, expr)
        return True, expr

    def beforeChild(self, node, child):
        child_type = child.__class__
        if child_type in native_types:
            return False, (self.CONSTANT, child)
        if not child.is_expression_type():
            if child.is_fixed():
                return False, (self.CONSTANT, value(child))
            else:
                return False, (self.VARIABLE, child)

        if child_type is not LinearExpression:
            return True, None

        # Because we are going to modify the LinearExpression in this
        # walker, we need to make a copy of the LinearExpression from
        # the orifinal expression tree.
        linear = self._get_linear()
        idMap = {}
        zeros = set()

        linear.constant = child.constant
        for c, v in zip(child.linear_coefs, child.linear_vars):
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
                    idMap[_id] = len(linear.linear_coefs)
                    linear.linear_vars.append(v)
                    linear.linear_coefs.append(c)
        if zeros:
            self._finalize_linear(zeros, linear)
        return False, (self.LINEAR, linear)

    def enterNode(self, node):
        if node._precedence() != SumExpressionBase.PRECEDENCE:
            return None, []
        return None, ({}, set(), self._get_linear())

    def acceptChildResult(self, node, data, child_result):
        if data.__class__ is list:
            # General expression... cache the child result until the end
            data.append(child_result)
        else:
            # Linear Expression
            idMap, zeros, linear = data
            child_type, child = child_result
            if not self._update_linear_expr(
                    idMap, zeros, linear, child_type, child):
                # General (nonlinear) sub-expression.  We will convert
                # this from a Linear expression to a general summation
                if linear.constant or linear.linear_vars:
                    if zeros:
                        self._finalize_linear(zeros, linear)
                    if not linear.linear_coefs:
                        const, linear.constant = linear.constant, 0
                        self.linearExprPool = (self.linearExprPool, linear)
                        return [(self.CONSTANT, const), child_result]
                    return [(self.LINEAR, linear), child_result]
                else:
                    return [child_result]

        return data

    def exitNode(self, node, data):
        if data.__class__ is tuple:
            # Linear expression
            _, zeros, linear = data
            if zeros:
                self._finalize_linear(zeros, linear)
            if not linear.linear_coefs:
                const, linear.constant = linear.constant, 0
                self.linearExprPool = (self.linearExprPool, linear)
                return (self.CONSTANT, const)
            return (self.LINEAR, linear)
        #
        # General (nonlinear) expression...
        #
        # If all the arguments are constants, then evaluate this node
        if all(_[0] is self.CONSTANT for _ in data):
            return (
                self.CONSTANT,
                node._apply_operation(tuple(_[1] for _ in data))
            )
        # We need special handling for Product/Division expressions
        if len(data) == 2:
            if isinstance(node, ProductExpression):
                if data[1][0] is self.CONSTANT:
                    data = (data[1], data[0])
                if data[0][0] is self.CONSTANT:
                    if not data[0][1]:
                        return (self.CONSTANT, 0)
                    elif data[1][0] is self.VARIABLE:
                        return (
                            self.MONOMIAL,
                            self._get_monomial((data[0][1], data[1][1]))
                        )
                    elif data[1][0] is self.MONOMIAL:
                        _args = data[1][1]._args_
                        data[1][1]._args_ = (_args[0] * data[0][1], _args[1])
                        return data[1]
                    elif data[1][0] is self.LINEAR:
                        mul = data[0][1]
                        data[1][1].constant *= mul
                        for i in xrange(len(data[1][1].linear_coefs)):
                            data[1][1].linear_coefs[i] *= mul
                        return data[1]
            elif isinstance(node, DivisionExpression):
                if data[1][0] is self.CONSTANT:
                    div = data[1][1]
                    if data[0][0] is self.VARIABLE:
                        return (
                            self.MONOMIAL,
                            self._get_monomial((1./div, data[0][1]))
                        )
                    elif data[0][0] is self.MONOMIAL:
                        _args = data[0][1]._args_
                        data[0][1]._args_ = (_args[0] / div, _args[1])
                        return data[0]
                    elif data[0][0] is self.LINEAR:
                        data[0][1].constant /= div
                        for i in xrange(len(data[0][1].linear_coefs)):
                            data[0][1].linear_coefs[i] /= div
                        return data[0]
        return (
            self.GENERAL,
            node.create_node_with_local_data(tuple(_[1] for _ in data))
        )

    def finalizeResult(self, result):
        result_type, expr = result
        ans = StandardRepn()
        if result_type is self.LINEAR:
            ans.constant = expr.constant
            ans.linear_vars = tuple(expr.linear_vars)
            ans.linear_coefs = tuple(expr.linear_coefs)
        elif result_type is self.GENERAL:
            if isinstance(expr, SumExpressionBase):
                linear_terms = []
                linear = self._get_linear()
                idMap = {}
                zeros = set()
                for i,term in enumerate(expr._args_):
                    term_type = type(term)
                    if term_type in native_types:
                        term_type = self.CONSTANT
                    elif not term.is_expression_type():
                        term_type = self.VARIABLE
                    elif term_type is MonomialTermExpression:
                        term_type = self.MONOMIAL
                    elif term_type is LinearExpression:
                        term_type = self.LINEAR
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
        elif result_type is self.MONOMIAL:
            c,v = expr._args_
            if c:
                ans.linear_coefs = (c,)
                ans.linear_vars = (v,)
            self.monomialExprPool = (self.monomialExprPool, expr)
        elif result_type is self.VARIABLE:
            ans.linear_coefs = (1,)
            ans.linear_vars = (expr,)
        else:
            ans.constant = expr

        return ans

    def _get_linear(self):
        if self.linearExprPool:
            self.linearExprPool, ans = self.linearExprPool
            return ans
        else:
            return LinearExpression()

    def _get_monomial(self, args):
        if self.monomialExprPool:
            self.monomialExprPool, ans = self.monomialExprPool
            ans._args_ = args
            return ans
        else:
            return MonomialTermExpression(args)

    def _finalize_linear(zeros, linear):
        for i in reversed(sorted(zeros)):
            if not linear.linear_coefs[i]:
                linear.linear_coefs.pop(i)
                linear.linear_vars.pop(i)

    def _update_linear_expr(
            self, idMap, zeros, linear, child_type, child):
        if child_type is self.MONOMIAL:
            self.monomialExprPool = (self.monomialExprPool, child)
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
        elif child_type is self.CONSTANT:
            linear.constant += child
        elif child_type is self.LINEAR:
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
