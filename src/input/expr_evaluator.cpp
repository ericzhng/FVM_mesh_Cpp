/**
 * @file expr_evaluator.cpp
 * @brief Implementation of the expression evaluator using recursive descent parsing.
 */

#include "expr_evaluator.hpp"

#include <cctype>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace fvm
{

    // =============================================================================
    // Token Types
    // =============================================================================

    enum class TokenType
    {
        Number,
        Variable, // X or Y
        Plus,
        Minus,
        Star,
        Slash,
        Caret, // Power
        LParen,
        RParen,
        Equal,    // ==
        NotEqual, // !=
        Less,
        Greater,
        LessEq,
        GreaterEq,
        And,  // &&
        Or,   // ||
        Not,  // !
        Func, // sqrt, abs, sin, cos
        End,
        Error
    };

    struct Token
    {
        TokenType type;
        std::string value;
        double numValue = 0.0;
    };

    // =============================================================================
    // Lexer
    // =============================================================================

    class Lexer
    {
    public:
        explicit Lexer(const std::string &input) : input_(input), pos_(0) {}

        Token nextToken()
        {
            skipWhitespace();

            if (pos_ >= input_.size())
            {
                return {TokenType::End, "", 0};
            }

            char c = input_[pos_];

            // Numbers
            if (std::isdigit(c) || (c == '.' && pos_ + 1 < input_.size() && std::isdigit(input_[pos_ + 1])))
            {
                return parseNumber();
            }

            // Identifiers (variables and functions)
            if (std::isalpha(c) || c == '_')
            {
                return parseIdentifier();
            }

            // Operators
            switch (c)
            {
            case '+':
                ++pos_;
                return {TokenType::Plus, "+", 0};
            case '-':
                ++pos_;
                return {TokenType::Minus, "-", 0};
            case '*':
                ++pos_;
                return {TokenType::Star, "*", 0};
            case '/':
                ++pos_;
                return {TokenType::Slash, "/", 0};
            case '^':
                ++pos_;
                return {TokenType::Caret, "^", 0};
            case '(':
                ++pos_;
                return {TokenType::LParen, "(", 0};
            case ')':
                ++pos_;
                return {TokenType::RParen, ")", 0};
            case '!':
                ++pos_;
                if (pos_ < input_.size() && input_[pos_] == '=')
                {
                    ++pos_;
                    return {TokenType::NotEqual, "!=", 0};
                }
                return {TokenType::Not, "!", 0};
            case '=':
                ++pos_;
                if (pos_ < input_.size() && input_[pos_] == '=')
                {
                    ++pos_;
                    return {TokenType::Equal, "==", 0};
                }
                return {TokenType::Error, "Single '=' not allowed, use '=='", 0};
            case '<':
                ++pos_;
                if (pos_ < input_.size() && input_[pos_] == '=')
                {
                    ++pos_;
                    return {TokenType::LessEq, "<=", 0};
                }
                return {TokenType::Less, "<", 0};
            case '>':
                ++pos_;
                if (pos_ < input_.size() && input_[pos_] == '=')
                {
                    ++pos_;
                    return {TokenType::GreaterEq, ">=", 0};
                }
                return {TokenType::Greater, ">", 0};
            case '&':
                ++pos_;
                if (pos_ < input_.size() && input_[pos_] == '&')
                {
                    ++pos_;
                    return {TokenType::And, "&&", 0};
                }
                return {TokenType::Error, "Single '&' not allowed, use '&&'", 0};
            case '|':
                ++pos_;
                if (pos_ < input_.size() && input_[pos_] == '|')
                {
                    ++pos_;
                    return {TokenType::Or, "||", 0};
                }
                return {TokenType::Error, "Single '|' not allowed, use '||'", 0};
            default:
                return {TokenType::Error, std::string("Unexpected character: ") + c, 0};
            }
        }

    private:
        void skipWhitespace()
        {
            while (pos_ < input_.size() && std::isspace(input_[pos_]))
            {
                ++pos_;
            }
        }

        Token parseNumber()
        {
            std::size_t start = pos_;
            while (pos_ < input_.size() && (std::isdigit(input_[pos_]) || input_[pos_] == '.'))
            {
                ++pos_;
            }
            // Handle scientific notation
            if (pos_ < input_.size() && (input_[pos_] == 'e' || input_[pos_] == 'E'))
            {
                ++pos_;
                if (pos_ < input_.size() && (input_[pos_] == '+' || input_[pos_] == '-'))
                {
                    ++pos_;
                }
                while (pos_ < input_.size() && std::isdigit(input_[pos_]))
                {
                    ++pos_;
                }
            }
            std::string numStr = input_.substr(start, pos_ - start);
            double value = std::stod(numStr);
            return {TokenType::Number, numStr, value};
        }

        Token parseIdentifier()
        {
            std::size_t start = pos_;
            while (pos_ < input_.size() && (std::isalnum(input_[pos_]) || input_[pos_] == '_'))
            {
                ++pos_;
            }
            std::string name = input_.substr(start, pos_ - start);

            // Check for variables
            if (name == "X" || name == "x")
            {
                return {TokenType::Variable, "X", 0};
            }
            if (name == "Y" || name == "y")
            {
                return {TokenType::Variable, "Y", 0};
            }

            // Check for functions
            if (name == "sqrt" || name == "abs" || name == "sin" || name == "cos" ||
                name == "tan" || name == "exp" || name == "log" || name == "pi")
            {
                return {TokenType::Func, name, 0};
            }

            return {TokenType::Error, "Unknown identifier: " + name, 0};
        }

        std::string input_;
        std::size_t pos_;
    };

    // =============================================================================
    // AST Nodes
    // =============================================================================

    class ASTNode
    {
    public:
        virtual ~ASTNode() = default;
        virtual double eval(double x, double y) const = 0;
    };

    using ASTPtr = std::unique_ptr<ASTNode>;

    class NumberNode : public ASTNode
    {
    public:
        explicit NumberNode(double value) : value_(value) {}
        double eval(double, double) const override { return value_; }

    private:
        double value_;
    };

    class VariableNode : public ASTNode
    {
    public:
        explicit VariableNode(char var) : var_(var) {}
        double eval(double x, double y) const override
        {
            return (var_ == 'X') ? x : y;
        }

    private:
        char var_;
    };

    class BinaryOpNode : public ASTNode
    {
    public:
        BinaryOpNode(TokenType op, ASTPtr left, ASTPtr right)
            : op_(op), left_(std::move(left)), right_(std::move(right)) {}

        double eval(double x, double y) const override
        {
            double lval = left_->eval(x, y);
            double rval = right_->eval(x, y);

            switch (op_)
            {
            case TokenType::Plus:
                return lval + rval;
            case TokenType::Minus:
                return lval - rval;
            case TokenType::Star:
                return lval * rval;
            case TokenType::Slash:
                return (rval != 0) ? lval / rval : 0;
            case TokenType::Caret:
                return std::pow(lval, rval);
            case TokenType::Equal:
                return (std::abs(lval - rval) < 1e-9) ? 1.0 : 0.0;
            case TokenType::NotEqual:
                return (std::abs(lval - rval) >= 1e-9) ? 1.0 : 0.0;
            case TokenType::Less:
                return (lval < rval) ? 1.0 : 0.0;
            case TokenType::Greater:
                return (lval > rval) ? 1.0 : 0.0;
            case TokenType::LessEq:
                return (lval <= rval) ? 1.0 : 0.0;
            case TokenType::GreaterEq:
                return (lval >= rval) ? 1.0 : 0.0;
            case TokenType::And:
                return (lval != 0 && rval != 0) ? 1.0 : 0.0;
            case TokenType::Or:
                return (lval != 0 || rval != 0) ? 1.0 : 0.0;
            default:
                return 0;
            }
        }

    private:
        TokenType op_;
        ASTPtr left_;
        ASTPtr right_;
    };

    class UnaryOpNode : public ASTNode
    {
    public:
        UnaryOpNode(TokenType op, ASTPtr operand)
            : op_(op), operand_(std::move(operand)) {}

        double eval(double x, double y) const override
        {
            double val = operand_->eval(x, y);
            switch (op_)
            {
            case TokenType::Minus:
                return -val;
            case TokenType::Not:
                return (val == 0) ? 1.0 : 0.0;
            default:
                return val;
            }
        }

    private:
        TokenType op_;
        ASTPtr operand_;
    };

    class FunctionNode : public ASTNode
    {
    public:
        FunctionNode(const std::string &func, ASTPtr arg)
            : func_(func), arg_(std::move(arg)) {}

        double eval(double x, double y) const override
        {
            double val = arg_ ? arg_->eval(x, y) : 0;

            if (func_ == "sqrt")
                return std::sqrt(val);
            if (func_ == "abs")
                return std::abs(val);
            if (func_ == "sin")
                return std::sin(val);
            if (func_ == "cos")
                return std::cos(val);
            if (func_ == "tan")
                return std::tan(val);
            if (func_ == "exp")
                return std::exp(val);
            if (func_ == "log")
                return std::log(val);
            if (func_ == "pi")
                return 3.14159265358979323846;

            return 0;
        }

    private:
        std::string func_;
        ASTPtr arg_;
    };

    // =============================================================================
    // Parser
    // =============================================================================

    class Parser
    {
    public:
        explicit Parser(const std::string &input) : lexer_(input)
        {
            advance();
        }

        ASTPtr parse()
        {
            auto result = parseOr();
            if (current_.type != TokenType::End)
            {
                error_ = "Unexpected token at end: " + current_.value;
                return nullptr;
            }
            return result;
        }

        const std::string &getError() const { return error_; }

    private:
        void advance()
        {
            current_ = lexer_.nextToken();
            if (current_.type == TokenType::Error)
            {
                error_ = current_.value;
            }
        }

        bool match(TokenType type)
        {
            if (current_.type == type)
            {
                advance();
                return true;
            }
            return false;
        }

        // Grammar (precedence low to high):
        // or     -> and ( '||' and )*
        // and    -> comp ( '&&' comp )*
        // comp   -> add ( ('==' | '!=' | '<' | '>' | '<=' | '>=') add )*
        // add    -> mult ( ('+' | '-') mult )*
        // mult   -> power ( ('*' | '/') power )*
        // power  -> unary ( '^' power )?
        // unary  -> ('-' | '!') unary | primary
        // primary -> number | variable | func '(' expr ')' | '(' expr ')'

        ASTPtr parseOr()
        {
            auto left = parseAnd();
            while (current_.type == TokenType::Or)
            {
                advance();
                auto right = parseAnd();
                left = std::make_unique<BinaryOpNode>(TokenType::Or, std::move(left), std::move(right));
            }
            return left;
        }

        ASTPtr parseAnd()
        {
            auto left = parseComparison();
            while (current_.type == TokenType::And)
            {
                advance();
                auto right = parseComparison();
                left = std::make_unique<BinaryOpNode>(TokenType::And, std::move(left), std::move(right));
            }
            return left;
        }

        ASTPtr parseComparison()
        {
            auto left = parseAdd();
            while (current_.type == TokenType::Equal ||
                   current_.type == TokenType::NotEqual ||
                   current_.type == TokenType::Less ||
                   current_.type == TokenType::Greater ||
                   current_.type == TokenType::LessEq ||
                   current_.type == TokenType::GreaterEq)
            {
                TokenType op = current_.type;
                advance();
                auto right = parseAdd();
                left = std::make_unique<BinaryOpNode>(op, std::move(left), std::move(right));
            }
            return left;
        }

        ASTPtr parseAdd()
        {
            auto left = parseMult();
            while (current_.type == TokenType::Plus || current_.type == TokenType::Minus)
            {
                TokenType op = current_.type;
                advance();
                auto right = parseMult();
                left = std::make_unique<BinaryOpNode>(op, std::move(left), std::move(right));
            }
            return left;
        }

        ASTPtr parseMult()
        {
            auto left = parsePower();
            while (current_.type == TokenType::Star || current_.type == TokenType::Slash)
            {
                TokenType op = current_.type;
                advance();
                auto right = parsePower();
                left = std::make_unique<BinaryOpNode>(op, std::move(left), std::move(right));
            }
            return left;
        }

        ASTPtr parsePower()
        {
            auto left = parseUnary();
            if (current_.type == TokenType::Caret)
            {
                advance();
                auto right = parsePower(); // Right-associative
                left = std::make_unique<BinaryOpNode>(TokenType::Caret, std::move(left), std::move(right));
            }
            return left;
        }

        ASTPtr parseUnary()
        {
            if (current_.type == TokenType::Minus)
            {
                advance();
                auto operand = parseUnary();
                return std::make_unique<UnaryOpNode>(TokenType::Minus, std::move(operand));
            }
            if (current_.type == TokenType::Not)
            {
                advance();
                auto operand = parseUnary();
                return std::make_unique<UnaryOpNode>(TokenType::Not, std::move(operand));
            }
            return parsePrimary();
        }

        ASTPtr parsePrimary()
        {
            if (current_.type == TokenType::Number)
            {
                double val = current_.numValue;
                advance();
                return std::make_unique<NumberNode>(val);
            }

            if (current_.type == TokenType::Variable)
            {
                char var = current_.value[0];
                advance();
                return std::make_unique<VariableNode>(var);
            }

            if (current_.type == TokenType::Func)
            {
                std::string func = current_.value;
                advance();

                // 'pi' is a constant, no arguments
                if (func == "pi")
                {
                    return std::make_unique<FunctionNode>(func, nullptr);
                }

                // Function call with argument
                if (!match(TokenType::LParen))
                {
                    error_ = "Expected '(' after function " + func;
                    return nullptr;
                }
                auto arg = parseOr();
                if (!match(TokenType::RParen))
                {
                    error_ = "Expected ')' after function argument";
                    return nullptr;
                }
                return std::make_unique<FunctionNode>(func, std::move(arg));
            }

            if (current_.type == TokenType::LParen)
            {
                advance();
                auto expr = parseOr();
                if (!match(TokenType::RParen))
                {
                    error_ = "Expected ')'";
                    return nullptr;
                }
                return expr;
            }

            error_ = "Unexpected token: " + current_.value;
            return nullptr;
        }

        Lexer lexer_;
        Token current_;
        std::string error_;
    };

    // =============================================================================
    // ExprEvaluator Implementation
    // =============================================================================

    struct ExprEvaluator::Impl
    {
        std::string expression;
        ASTPtr ast;
        bool valid = false;
        std::string error;
    };

    ExprEvaluator::ExprEvaluator(const std::string &expression)
        : pImpl_(std::make_unique<Impl>())
    {
        pImpl_->expression = expression;

        Parser parser(expression);
        pImpl_->ast = parser.parse();

        if (pImpl_->ast && parser.getError().empty())
        {
            pImpl_->valid = true;
        }
        else
        {
            pImpl_->valid = false;
            pImpl_->error = parser.getError();
            if (pImpl_->error.empty())
            {
                pImpl_->error = "Failed to parse expression";
            }
        }
    }

    ExprEvaluator::~ExprEvaluator() = default;

    ExprEvaluator::ExprEvaluator(ExprEvaluator &&other) noexcept = default;

    ExprEvaluator &ExprEvaluator::operator=(ExprEvaluator &&other) noexcept = default;

    double ExprEvaluator::evaluate(double x, double y) const
    {
        if (!pImpl_->valid || !pImpl_->ast)
        {
            return 0.0;
        }
        return pImpl_->ast->eval(x, y);
    }

    bool ExprEvaluator::matches(double x, double y) const
    {
        return evaluate(x, y) != 0.0;
    }

    bool ExprEvaluator::isValid() const
    {
        return pImpl_->valid;
    }

    const std::string &ExprEvaluator::getError() const
    {
        return pImpl_->error;
    }

    const std::string &ExprEvaluator::getExpression() const
    {
        return pImpl_->expression;
    }

} // namespace fvm
