#pragma once

/**
 * @file expr_evaluator.hpp
 * @brief Simple expression evaluator for boundary condition expressions.
 *
 * This evaluator supports mathematical expressions with variables X and Y,
 * used to identify boundary edges based on their midpoint coordinates.
 *
 * Supported features:
 *   - Variables: X, Y
 *   - Numbers: integers and floating-point
 *   - Arithmetic: +, -, *, /, ^ (power)
 *   - Comparison: ==, !=, <, >, <=, >=
 *   - Logical: ||, &&, !
 *   - Functions: sqrt(), abs(), sin(), cos()
 *   - Parentheses for grouping
 *
 * Examples:
 *   - "X == 0"
 *   - "Y >= 0 && Y <= 1"
 *   - "X^2 + Y^2 <= 0.25"
 *   - "abs(X - 0.5) < 0.01"
 */

#include "common/fvm_export.hpp"

#include <memory>
#include <string>

namespace fvm
{

    /**
     * @brief Expression evaluator for 2D coordinate-based expressions.
     *
     * The evaluator parses an expression string once during construction,
     * then can efficiently evaluate it for different (X, Y) values.
     */
    class FVM_API ExprEvaluator
    {
    public:
        /**
         * @brief Construct an evaluator for the given expression.
         * @param expression The mathematical expression string
         */
        explicit ExprEvaluator(const std::string &expression);

        /// Destructor
        ~ExprEvaluator();

        /// Move constructor
        ExprEvaluator(ExprEvaluator &&other) noexcept;

        /// Move assignment
        ExprEvaluator &operator=(ExprEvaluator &&other) noexcept;

        // Disable copy (expression tree is not trivially copyable)
        ExprEvaluator(const ExprEvaluator &) = delete;
        ExprEvaluator &operator=(const ExprEvaluator &) = delete;

        /**
         * @brief Evaluate the expression at given coordinates.
         * @param x X-coordinate value
         * @param y Y-coordinate value
         * @return Result of expression evaluation (non-zero = true for comparisons)
         */
        double evaluate(double x, double y) const;

        /**
         * @brief Evaluate expression as boolean (for boundary matching).
         * @param x X-coordinate value
         * @param y Y-coordinate value
         * @return true if expression evaluates to non-zero
         */
        bool matches(double x, double y) const;

        /**
         * @brief Check if the expression was parsed successfully.
         * @return true if expression is valid
         */
        bool isValid() const;

        /**
         * @brief Get the error message if parsing failed.
         * @return Error message, or empty string if valid
         */
        const std::string &getError() const;

        /**
         * @brief Get the original expression string.
         * @return The expression that was parsed
         */
        const std::string &getExpression() const;

    private:
        struct Impl;
        std::unique_ptr<Impl> pImpl_;
    };

} // namespace fvm
