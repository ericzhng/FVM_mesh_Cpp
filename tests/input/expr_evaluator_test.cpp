#include <gtest/gtest.h>
#include "input/expr_evaluator.hpp"
#include <cmath>

namespace fvm
{
    namespace testing
    {

        // =============================================================================
        // Basic Arithmetic Tests
        // =============================================================================

        TEST(ExprEvaluatorArithmetic, NumberLiteral)
        {
            ExprEvaluator eval("42");
            EXPECT_TRUE(eval.isValid());
            EXPECT_DOUBLE_EQ(eval.evaluate(0, 0), 42.0);
        }

        TEST(ExprEvaluatorArithmetic, FloatLiteral)
        {
            ExprEvaluator eval("3.14159");
            EXPECT_TRUE(eval.isValid());
            EXPECT_NEAR(eval.evaluate(0, 0), 3.14159, 1e-10);
        }

        TEST(ExprEvaluatorArithmetic, ScientificNotation)
        {
            ExprEvaluator eval("1.5e-3");
            EXPECT_TRUE(eval.isValid());
            EXPECT_NEAR(eval.evaluate(0, 0), 0.0015, 1e-15);
        }

        TEST(ExprEvaluatorArithmetic, Addition)
        {
            ExprEvaluator eval("3 + 5");
            EXPECT_TRUE(eval.isValid());
            EXPECT_DOUBLE_EQ(eval.evaluate(0, 0), 8.0);
        }

        TEST(ExprEvaluatorArithmetic, Subtraction)
        {
            ExprEvaluator eval("10 - 4");
            EXPECT_TRUE(eval.isValid());
            EXPECT_DOUBLE_EQ(eval.evaluate(0, 0), 6.0);
        }

        TEST(ExprEvaluatorArithmetic, Multiplication)
        {
            ExprEvaluator eval("6 * 7");
            EXPECT_TRUE(eval.isValid());
            EXPECT_DOUBLE_EQ(eval.evaluate(0, 0), 42.0);
        }

        TEST(ExprEvaluatorArithmetic, Division)
        {
            ExprEvaluator eval("15 / 3");
            EXPECT_TRUE(eval.isValid());
            EXPECT_DOUBLE_EQ(eval.evaluate(0, 0), 5.0);
        }

        TEST(ExprEvaluatorArithmetic, Power)
        {
            ExprEvaluator eval("2 ^ 3");
            EXPECT_TRUE(eval.isValid());
            EXPECT_DOUBLE_EQ(eval.evaluate(0, 0), 8.0);
        }

        TEST(ExprEvaluatorArithmetic, UnaryMinus)
        {
            ExprEvaluator eval("-5");
            EXPECT_TRUE(eval.isValid());
            EXPECT_DOUBLE_EQ(eval.evaluate(0, 0), -5.0);
        }

        TEST(ExprEvaluatorArithmetic, OperatorPrecedence)
        {
            ExprEvaluator eval("2 + 3 * 4");
            EXPECT_TRUE(eval.isValid());
            EXPECT_DOUBLE_EQ(eval.evaluate(0, 0), 14.0); // 2 + (3*4) = 14
        }

        TEST(ExprEvaluatorArithmetic, Parentheses)
        {
            ExprEvaluator eval("(2 + 3) * 4");
            EXPECT_TRUE(eval.isValid());
            EXPECT_DOUBLE_EQ(eval.evaluate(0, 0), 20.0);
        }

        TEST(ExprEvaluatorArithmetic, NestedParentheses)
        {
            ExprEvaluator eval("((1 + 2) * (3 + 4))");
            EXPECT_TRUE(eval.isValid());
            EXPECT_DOUBLE_EQ(eval.evaluate(0, 0), 21.0);
        }

        TEST(ExprEvaluatorArithmetic, PowerRightAssociative)
        {
            ExprEvaluator eval("2 ^ 3 ^ 2");
            EXPECT_TRUE(eval.isValid());
            EXPECT_DOUBLE_EQ(eval.evaluate(0, 0), 512.0); // 2^(3^2) = 2^9 = 512
        }

        // =============================================================================
        // Variable Tests
        // =============================================================================

        TEST(ExprEvaluatorVariables, VariableX)
        {
            ExprEvaluator eval("X");
            EXPECT_TRUE(eval.isValid());
            EXPECT_DOUBLE_EQ(eval.evaluate(5.0, 3.0), 5.0);
        }

        TEST(ExprEvaluatorVariables, VariableY)
        {
            ExprEvaluator eval("Y");
            EXPECT_TRUE(eval.isValid());
            EXPECT_DOUBLE_EQ(eval.evaluate(5.0, 3.0), 3.0);
        }

        TEST(ExprEvaluatorVariables, LowercaseX)
        {
            ExprEvaluator eval("x");
            EXPECT_TRUE(eval.isValid());
            EXPECT_DOUBLE_EQ(eval.evaluate(7.0, 2.0), 7.0);
        }

        TEST(ExprEvaluatorVariables, LowercaseY)
        {
            ExprEvaluator eval("y");
            EXPECT_TRUE(eval.isValid());
            EXPECT_DOUBLE_EQ(eval.evaluate(7.0, 2.0), 2.0);
        }

        TEST(ExprEvaluatorVariables, ExpressionWithVariables)
        {
            ExprEvaluator eval("X + Y");
            EXPECT_TRUE(eval.isValid());
            EXPECT_DOUBLE_EQ(eval.evaluate(3.0, 4.0), 7.0);
        }

        TEST(ExprEvaluatorVariables, ComplexExpression)
        {
            ExprEvaluator eval("X^2 + Y^2");
            EXPECT_TRUE(eval.isValid());
            EXPECT_DOUBLE_EQ(eval.evaluate(3.0, 4.0), 25.0);
        }

        // =============================================================================
        // Comparison Tests
        // =============================================================================

        TEST(ExprEvaluatorComparison, Equal)
        {
            ExprEvaluator eval("X == 0");
            EXPECT_TRUE(eval.isValid());
            EXPECT_TRUE(eval.matches(0.0, 0.0));
            EXPECT_FALSE(eval.matches(1.0, 0.0));
        }

        TEST(ExprEvaluatorComparison, NotEqual)
        {
            ExprEvaluator eval("X != 0");
            EXPECT_TRUE(eval.isValid());
            EXPECT_FALSE(eval.matches(0.0, 0.0));
            EXPECT_TRUE(eval.matches(1.0, 0.0));
        }

        TEST(ExprEvaluatorComparison, LessThan)
        {
            ExprEvaluator eval("X < 5");
            EXPECT_TRUE(eval.isValid());
            EXPECT_TRUE(eval.matches(3.0, 0.0));
            EXPECT_FALSE(eval.matches(5.0, 0.0));
            EXPECT_FALSE(eval.matches(7.0, 0.0));
        }

        TEST(ExprEvaluatorComparison, GreaterThan)
        {
            ExprEvaluator eval("X > 5");
            EXPECT_TRUE(eval.isValid());
            EXPECT_FALSE(eval.matches(3.0, 0.0));
            EXPECT_FALSE(eval.matches(5.0, 0.0));
            EXPECT_TRUE(eval.matches(7.0, 0.0));
        }

        TEST(ExprEvaluatorComparison, LessOrEqual)
        {
            ExprEvaluator eval("X <= 5");
            EXPECT_TRUE(eval.isValid());
            EXPECT_TRUE(eval.matches(3.0, 0.0));
            EXPECT_TRUE(eval.matches(5.0, 0.0));
            EXPECT_FALSE(eval.matches(7.0, 0.0));
        }

        TEST(ExprEvaluatorComparison, GreaterOrEqual)
        {
            ExprEvaluator eval("X >= 5");
            EXPECT_TRUE(eval.isValid());
            EXPECT_FALSE(eval.matches(3.0, 0.0));
            EXPECT_TRUE(eval.matches(5.0, 0.0));
            EXPECT_TRUE(eval.matches(7.0, 0.0));
        }

        // =============================================================================
        // Logical Operator Tests
        // =============================================================================

        TEST(ExprEvaluatorLogical, And)
        {
            ExprEvaluator eval("X >= 0 && X <= 1");
            EXPECT_TRUE(eval.isValid());
            EXPECT_FALSE(eval.matches(-1.0, 0.0));
            EXPECT_TRUE(eval.matches(0.0, 0.0));
            EXPECT_TRUE(eval.matches(0.5, 0.0));
            EXPECT_TRUE(eval.matches(1.0, 0.0));
            EXPECT_FALSE(eval.matches(2.0, 0.0));
        }

        TEST(ExprEvaluatorLogical, Or)
        {
            ExprEvaluator eval("X == 0 || Y == 0");
            EXPECT_TRUE(eval.isValid());
            EXPECT_TRUE(eval.matches(0.0, 5.0));
            EXPECT_TRUE(eval.matches(5.0, 0.0));
            EXPECT_TRUE(eval.matches(0.0, 0.0));
            EXPECT_FALSE(eval.matches(1.0, 1.0));
        }

        TEST(ExprEvaluatorLogical, Not)
        {
            ExprEvaluator eval("!(X == 0)");
            EXPECT_TRUE(eval.isValid());
            EXPECT_FALSE(eval.matches(0.0, 0.0));
            EXPECT_TRUE(eval.matches(1.0, 0.0));
        }

        TEST(ExprEvaluatorLogical, ComplexLogical)
        {
            ExprEvaluator eval("(X >= 0 && X <= 1) && (Y >= 0 && Y <= 1)");
            EXPECT_TRUE(eval.isValid());
            EXPECT_TRUE(eval.matches(0.5, 0.5));
            EXPECT_FALSE(eval.matches(0.5, 1.5));
            EXPECT_FALSE(eval.matches(-0.5, 0.5));
        }

        // =============================================================================
        // Function Tests
        // =============================================================================

        TEST(ExprEvaluatorFunctions, Sqrt)
        {
            ExprEvaluator eval("sqrt(16)");
            EXPECT_TRUE(eval.isValid());
            EXPECT_DOUBLE_EQ(eval.evaluate(0, 0), 4.0);
        }

        TEST(ExprEvaluatorFunctions, Abs)
        {
            ExprEvaluator eval("abs(-5)");
            EXPECT_TRUE(eval.isValid());
            EXPECT_DOUBLE_EQ(eval.evaluate(0, 0), 5.0);
        }

        TEST(ExprEvaluatorFunctions, AbsWithVariable)
        {
            ExprEvaluator eval("abs(X - 0.5)");
            EXPECT_TRUE(eval.isValid());
            EXPECT_DOUBLE_EQ(eval.evaluate(0.3, 0), 0.2);
            EXPECT_DOUBLE_EQ(eval.evaluate(0.7, 0), 0.2);
        }

        TEST(ExprEvaluatorFunctions, Sin)
        {
            ExprEvaluator eval("sin(0)");
            EXPECT_TRUE(eval.isValid());
            EXPECT_NEAR(eval.evaluate(0, 0), 0.0, 1e-10);
        }

        TEST(ExprEvaluatorFunctions, Cos)
        {
            ExprEvaluator eval("cos(0)");
            EXPECT_TRUE(eval.isValid());
            EXPECT_NEAR(eval.evaluate(0, 0), 1.0, 1e-10);
        }

        TEST(ExprEvaluatorFunctions, Tan)
        {
            ExprEvaluator eval("tan(0)");
            EXPECT_TRUE(eval.isValid());
            EXPECT_NEAR(eval.evaluate(0, 0), 0.0, 1e-10);
        }

        TEST(ExprEvaluatorFunctions, Exp)
        {
            ExprEvaluator eval("exp(0)");
            EXPECT_TRUE(eval.isValid());
            EXPECT_DOUBLE_EQ(eval.evaluate(0, 0), 1.0);
        }

        TEST(ExprEvaluatorFunctions, Log)
        {
            ExprEvaluator eval("log(1)");
            EXPECT_TRUE(eval.isValid());
            EXPECT_DOUBLE_EQ(eval.evaluate(0, 0), 0.0);
        }

        TEST(ExprEvaluatorFunctions, Pi)
        {
            ExprEvaluator eval("pi");
            EXPECT_TRUE(eval.isValid());
            EXPECT_NEAR(eval.evaluate(0, 0), 3.14159265358979323846, 1e-14);
        }

        TEST(ExprEvaluatorFunctions, SinWithPi)
        {
            ExprEvaluator eval("sin(pi)");
            EXPECT_TRUE(eval.isValid());
            EXPECT_NEAR(eval.evaluate(0, 0), 0.0, 1e-10);
        }

        TEST(ExprEvaluatorFunctions, NestedFunction)
        {
            ExprEvaluator eval("sqrt(abs(-16))");
            EXPECT_TRUE(eval.isValid());
            EXPECT_DOUBLE_EQ(eval.evaluate(0, 0), 4.0);
        }

        // =============================================================================
        // Boundary Condition Expression Tests
        // =============================================================================

        TEST(ExprEvaluatorBoundary, LeftBoundary)
        {
            ExprEvaluator eval("X == 0");
            EXPECT_TRUE(eval.isValid());
            EXPECT_TRUE(eval.matches(0.0, 0.5));
            EXPECT_FALSE(eval.matches(1.0, 0.5));
        }

        TEST(ExprEvaluatorBoundary, RightBoundary)
        {
            ExprEvaluator eval("X == 1");
            EXPECT_TRUE(eval.isValid());
            EXPECT_TRUE(eval.matches(1.0, 0.5));
            EXPECT_FALSE(eval.matches(0.0, 0.5));
        }

        TEST(ExprEvaluatorBoundary, BottomBoundary)
        {
            ExprEvaluator eval("Y == 0");
            EXPECT_TRUE(eval.isValid());
            EXPECT_TRUE(eval.matches(0.5, 0.0));
            EXPECT_FALSE(eval.matches(0.5, 1.0));
        }

        TEST(ExprEvaluatorBoundary, TopBoundary)
        {
            ExprEvaluator eval("Y == 1");
            EXPECT_TRUE(eval.isValid());
            EXPECT_TRUE(eval.matches(0.5, 1.0));
            EXPECT_FALSE(eval.matches(0.5, 0.0));
        }

        TEST(ExprEvaluatorBoundary, CircularBoundary)
        {
            ExprEvaluator eval("X^2 + Y^2 <= 0.25");
            EXPECT_TRUE(eval.isValid());
            EXPECT_TRUE(eval.matches(0.0, 0.0));   // Center
            EXPECT_TRUE(eval.matches(0.25, 0.25)); // Inside (0.125 <= 0.25)
            EXPECT_FALSE(eval.matches(0.5, 0.5));  // Outside (0.5 > 0.25)
        }

        TEST(ExprEvaluatorBoundary, RangeCheck)
        {
            ExprEvaluator eval("Y >= 0 && Y <= 1");
            EXPECT_TRUE(eval.isValid());
            EXPECT_TRUE(eval.matches(5.0, 0.5));
            EXPECT_FALSE(eval.matches(5.0, -0.1));
            EXPECT_FALSE(eval.matches(5.0, 1.1));
        }

        TEST(ExprEvaluatorBoundary, ToleranceCheck)
        {
            ExprEvaluator eval("abs(X - 0.5) < 0.01");
            EXPECT_TRUE(eval.isValid());
            EXPECT_TRUE(eval.matches(0.5, 0.0));
            EXPECT_TRUE(eval.matches(0.505, 0.0));
            EXPECT_FALSE(eval.matches(0.52, 0.0));
        }

        // =============================================================================
        // Error Handling Tests
        // =============================================================================

        TEST(ExprEvaluatorErrors, InvalidExpression)
        {
            ExprEvaluator eval("X + ");
            EXPECT_FALSE(eval.isValid());
            EXPECT_FALSE(eval.getError().empty());
        }

        TEST(ExprEvaluatorErrors, UnknownIdentifier)
        {
            ExprEvaluator eval("Z + 1");
            EXPECT_FALSE(eval.isValid());
            EXPECT_FALSE(eval.getError().empty());
        }

        TEST(ExprEvaluatorErrors, SingleEqualSign)
        {
            ExprEvaluator eval("X = 0");
            EXPECT_FALSE(eval.isValid());
            EXPECT_FALSE(eval.getError().empty());
        }

        TEST(ExprEvaluatorErrors, SingleAmpersand)
        {
            ExprEvaluator eval("X > 0 & Y > 0");
            EXPECT_FALSE(eval.isValid());
            EXPECT_FALSE(eval.getError().empty());
        }

        TEST(ExprEvaluatorErrors, SinglePipe)
        {
            ExprEvaluator eval("X == 0 | Y == 0");
            EXPECT_FALSE(eval.isValid());
            EXPECT_FALSE(eval.getError().empty());
        }

        TEST(ExprEvaluatorErrors, MismatchedParentheses)
        {
            ExprEvaluator eval("(X + Y");
            EXPECT_FALSE(eval.isValid());
            EXPECT_FALSE(eval.getError().empty());
        }

        TEST(ExprEvaluatorErrors, EmptyExpression)
        {
            ExprEvaluator eval("");
            EXPECT_FALSE(eval.isValid());
        }

        TEST(ExprEvaluatorErrors, InvalidReturnZeroIfInvalid)
        {
            ExprEvaluator eval("invalid++expression");
            EXPECT_FALSE(eval.isValid());
            EXPECT_DOUBLE_EQ(eval.evaluate(1.0, 1.0), 0.0);
        }

        // =============================================================================
        // Move Semantics Tests
        // =============================================================================

        TEST(ExprEvaluatorMoveSemantics, MoveConstruct)
        {
            ExprEvaluator eval1("X + Y");
            EXPECT_TRUE(eval1.isValid());

            ExprEvaluator eval2(std::move(eval1));
            EXPECT_TRUE(eval2.isValid());
            EXPECT_DOUBLE_EQ(eval2.evaluate(3.0, 4.0), 7.0);
        }

        TEST(ExprEvaluatorMoveSemantics, MoveAssign)
        {
            ExprEvaluator eval1("X + Y");
            ExprEvaluator eval2("X * Y");

            eval2 = std::move(eval1);
            EXPECT_TRUE(eval2.isValid());
            EXPECT_DOUBLE_EQ(eval2.evaluate(3.0, 4.0), 7.0);
        }

        // =============================================================================
        // GetExpression Tests
        // =============================================================================

        TEST(ExprEvaluatorMisc, GetExpression)
        {
            ExprEvaluator eval("X^2 + Y^2");
            EXPECT_EQ(eval.getExpression(), "X^2 + Y^2");
        }

        TEST(ExprEvaluatorMisc, WhitespaceHandling)
        {
            ExprEvaluator eval("   X   +   Y   ");
            EXPECT_TRUE(eval.isValid());
            EXPECT_DOUBLE_EQ(eval.evaluate(1.0, 2.0), 3.0);
        }

    } // namespace testing
} // namespace fvm
