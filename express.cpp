#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <memory>
#include <sstream>
#include <iomanip>
#include <cctype>
#include <algorithm>
#include <stdexcept>
#include <stack>
#include <functional>
#include <queue>

// Forward declarations
class ExpressionNode;
class Evaluator;

// Map for variable storage
std::map<std::string, double> variables;

// Function type definition
using MathFunction = std::function<double(const std::vector<double>&)>;

// Function registry
std::map<std::string, std::pair<MathFunction, int>> functions;

// Initialize built-in functions
void initFunctions() {
    functions["sin"] = {[](const std::vector<double>& args) { return std::sin(args[0]); }, 1};
    functions["cos"] = {[](const std::vector<double>& args) { return std::cos(args[0]); }, 1};
    functions["tan"] = {[](const std::vector<double>& args) { return std::tan(args[0]); }, 1};
    functions["sqrt"] = {[](const std::vector<double>& args) { return std::sqrt(args[0]); }, 1};
    functions["log"] = {[](const std::vector<double>& args) { return std::log(args[0]); }, 1};
    functions["exp"] = {[](const std::vector<double>& args) { return std::exp(args[0]); }, 1};
    functions["pow"] = {[](const std::vector<double>& args) { return std::pow(args[0], args[1]); }, 2};
    functions["max"] = {[](const std::vector<double>& args) { return std::max(args[0], args[1]); }, 2};
    functions["min"] = {[](const std::vector<double>& args) { return std::min(args[0], args[1]); }, 2};
    functions["abs"] = {[](const std::vector<double>& args) { return std::abs(args[0]); }, 1};
}

// Node type enumeration
enum class NodeType {
    CONSTANT,
    VARIABLE,
    UNARY_OPERATOR,
    BINARY_OPERATOR,
    FUNCTION
};

// Base Expression Node class
class ExpressionNode {
public:
    virtual ~ExpressionNode() = default;
    virtual double evaluate() const = 0;
    virtual std::string toString() const = 0;
    virtual NodeType getType() const = 0;
    virtual std::shared_ptr<ExpressionNode> clone() const = 0;
    virtual std::shared_ptr<ExpressionNode> simplify() = 0;
    virtual bool isConstant() const { return false; }
    virtual double getConstantValue() const { throw std::runtime_error("Not a constant node"); }
};

// Constant node implementation
class ConstantNode : public ExpressionNode {
private:
    double value;

public:
    ConstantNode(double val) : value(val) {}

    double evaluate() const override {
        return value;
    }

    std::string toString() const override {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(4);
        oss << value;
        std::string result = oss.str();
        // Remove trailing zeros
        if (result.find('.') != std::string::npos) {
            result = result.substr(0, result.find_last_not_of('0') + 1);
            if (result.back() == '.') {
                result.pop_back();
            }
        }
        return result;
    }

    NodeType getType() const override {
        return NodeType::CONSTANT;
    }

    std::shared_ptr<ExpressionNode> clone() const override {
        return std::make_shared<ConstantNode>(value);
    }

    std::shared_ptr<ExpressionNode> simplify() override {
        return clone();
    }

    bool isConstant() const override { return true; }
    double getConstantValue() const override { return value; }
};

// Variable node implementation
class VariableNode : public ExpressionNode {
private:
    std::string name;

public:
    VariableNode(const std::string& varName) : name(varName) {}

    double evaluate() const override {
        auto it = variables.find(name);
        if (it != variables.end()) {
            return it->second;
        }
        throw std::runtime_error("Variable '" + name + "' not defined");
    }

    std::string toString() const override {
        return name;
    }

    NodeType getType() const override {
        return NodeType::VARIABLE;
    }

    std::shared_ptr<ExpressionNode> clone() const override {
        return std::make_shared<VariableNode>(name);
    }

    std::shared_ptr<ExpressionNode> simplify() override {
        return clone();
    }

    const std::string& getName() const {
        return name;
    }
};

// Unary operator node implementation
class UnaryOperatorNode : public ExpressionNode {
private:
    char op;
    std::shared_ptr<ExpressionNode> operand;

public:
    UnaryOperatorNode(char oper, std::shared_ptr<ExpressionNode> expr)
        : op(oper), operand(std::move(expr)) {}

    double evaluate() const override {
        double val = operand->evaluate();
        switch (op) {
            case '-': return -val;
            case '+': return val;
            default: throw std::runtime_error("Unknown unary operator");
        }
    }

    std::string toString() const override {
        return op + std::string("(") + operand->toString() + ")";
    }

    NodeType getType() const override {
        return NodeType::UNARY_OPERATOR;
    }

    std::shared_ptr<ExpressionNode> clone() const override {
        return std::make_shared<UnaryOperatorNode>(op, operand->clone());
    }

    std::shared_ptr<ExpressionNode> simplify() override {
        auto simplifiedOperand = operand->simplify();

        // If operand is a constant, compute the result
        if (simplifiedOperand->isConstant()) {
            double val = simplifiedOperand->getConstantValue();
            switch (op) {
                case '-': return std::make_shared<ConstantNode>(-val);
                case '+': return simplifiedOperand; // +x = x
                default: throw std::runtime_error("Unknown unary operator");
            }
        }

        // For unary plus, just return the operand
        if (op == '+') {
            return simplifiedOperand;
        }

        // For unary minus of unary minus, remove double negation: -(-x) = x
        if (op == '-' && 
            simplifiedOperand->getType() == NodeType::UNARY_OPERATOR) {
            auto unaryNode = std::dynamic_pointer_cast<UnaryOperatorNode>(simplifiedOperand);
            if (unaryNode && unaryNode->getOperator() == '-') {
                return unaryNode->getOperand()->clone();
            }
        }

        return std::make_shared<UnaryOperatorNode>(op, simplifiedOperand);
    }

    char getOperator() const {
        return op;
    }

    std::shared_ptr<ExpressionNode> getOperand() const {
        return operand;
    }
};

// Binary operator node implementation
class BinaryOperatorNode : public ExpressionNode {
private:
    char op;
    std::shared_ptr<ExpressionNode> left;
    std::shared_ptr<ExpressionNode> right;

public:
    BinaryOperatorNode(char oper, std::shared_ptr<ExpressionNode> lhs, std::shared_ptr<ExpressionNode> rhs)
        : op(oper), left(std::move(lhs)), right(std::move(rhs)) {}

    double evaluate() const override {
        double lval = left->evaluate();
        double rval = right->evaluate();

        switch (op) {
            case '+': return lval + rval;
            case '-': return lval - rval;
            case '*': return lval * rval;
            case '/': 
                if (rval == 0) {
                    throw std::runtime_error("Division by zero");
                }
                return lval / rval;
            case '^': return std::pow(lval, rval);
            default: throw std::runtime_error("Unknown binary operator");
        }
    }

    std::string toString() const override {
        return "(" + left->toString() + " " + op + " " + right->toString() + ")";
    }

    NodeType getType() const override {
        return NodeType::BINARY_OPERATOR;
    }

    std::shared_ptr<ExpressionNode> clone() const override {
        return std::make_shared<BinaryOperatorNode>(op, left->clone(), right->clone());
    }

    std::shared_ptr<ExpressionNode> simplify() override {
        auto simplifiedLeft = left->simplify();
        auto simplifiedRight = right->simplify();

        // Both operands are constants, compute the result
        if (simplifiedLeft->isConstant() && simplifiedRight->isConstant()) {
            double lval = simplifiedLeft->getConstantValue();
            double rval = simplifiedRight->getConstantValue();

            switch (op) {
                case '+': return std::make_shared<ConstantNode>(lval + rval);
                case '-': return std::make_shared<ConstantNode>(lval - rval);
                case '*': return std::make_shared<ConstantNode>(lval * rval);
                case '/': 
                    if (rval == 0) {
                        throw std::runtime_error("Division by zero during simplification");
                    }
                    return std::make_shared<ConstantNode>(lval / rval);
                case '^': return std::make_shared<ConstantNode>(std::pow(lval, rval));
                default: throw std::runtime_error("Unknown binary operator");
            }
        }

        // Identify common algebraic simplifications
        switch (op) {
            case '+': {
                // 0 + x = x
                if (simplifiedLeft->isConstant() && simplifiedLeft->getConstantValue() == 0) {
                    return simplifiedRight;
                }
                // x + 0 = x
                if (simplifiedRight->isConstant() && simplifiedRight->getConstantValue() == 0) {
                    return simplifiedLeft;
                }
                break;
            }
            case '-': {
                // x - 0 = x
                if (simplifiedRight->isConstant() && simplifiedRight->getConstantValue() == 0) {
                    return simplifiedLeft;
                }
                // x - x = 0
                if (simplifiedLeft->toString() == simplifiedRight->toString()) {
                    return std::make_shared<ConstantNode>(0);
                }
                break;
            }
            case '*': {
                // 0 * x = 0
                if ((simplifiedLeft->isConstant() && simplifiedLeft->getConstantValue() == 0) ||
                    (simplifiedRight->isConstant() && simplifiedRight->getConstantValue() == 0)) {
                    return std::make_shared<ConstantNode>(0);
                }
                // 1 * x = x
                if (simplifiedLeft->isConstant() && simplifiedLeft->getConstantValue() == 1) {
                    return simplifiedRight;
                }
                // x * 1 = x
                if (simplifiedRight->isConstant() && simplifiedRight->getConstantValue() == 1) {
                    return simplifiedLeft;
                }
                break;
            }
            case '/': {
                // 0 / x = 0 (if x != 0)
                if (simplifiedLeft->isConstant() && simplifiedLeft->getConstantValue() == 0 &&
                    !(simplifiedRight->isConstant() && simplifiedRight->getConstantValue() == 0)) {
                    return std::make_shared<ConstantNode>(0);
                }
                // x / 1 = x
                if (simplifiedRight->isConstant() && simplifiedRight->getConstantValue() == 1) {
                    return simplifiedLeft;
                }
                // x / x = 1 (if x != 0)
                if (simplifiedLeft->toString() == simplifiedRight->toString()) {
                    return std::make_shared<ConstantNode>(1);
                }
                break;
            }
            case '^': {
                // x^0 = 1
                if (simplifiedRight->isConstant() && simplifiedRight->getConstantValue() == 0) {
                    return std::make_shared<ConstantNode>(1);
                }
                // x^1 = x
                if (simplifiedRight->isConstant() && simplifiedRight->getConstantValue() == 1) {
                    return simplifiedLeft;
                }
                // 1^x = 1
                if (simplifiedLeft->isConstant() && simplifiedLeft->getConstantValue() == 1) {
                    return std::make_shared<ConstantNode>(1);
                }
                break;
            }
        }

        return std::make_shared<BinaryOperatorNode>(op, simplifiedLeft, simplifiedRight);
    }

    char getOperator() const {
        return op;
    }

    std::shared_ptr<ExpressionNode> getLeft() const {
        return left;
    }

    std::shared_ptr<ExpressionNode> getRight() const {
        return right;
    }
};

// Function node implementation
class FunctionNode : public ExpressionNode {
private:
    std::string name;
    std::vector<std::shared_ptr<ExpressionNode>> arguments;

public:
    FunctionNode(const std::string& funcName, std::vector<std::shared_ptr<ExpressionNode>> args)
        : name(funcName), arguments(std::move(args)) {}

    double evaluate() const override {
        auto it = functions.find(name);
        if (it != functions.end()) {
            const auto& [func, expectedArgCount] = it->second;
            if (static_cast<int>(arguments.size()) != expectedArgCount) {
                throw std::runtime_error("Function '" + name + "' expects " +
                                         std::to_string(expectedArgCount) + " arguments, got " +
                                         std::to_string(arguments.size()));
            }

            std::vector<double> args;
            args.reserve(arguments.size());
            for (const auto& arg : arguments) {
                args.push_back(arg->evaluate());
            }
            return func(args);
        }
        throw std::runtime_error("Function '" + name + "' not defined");
    }

    std::string toString() const override {
        std::string result = name + "(";
        for (size_t i = 0; i < arguments.size(); ++i) {
            if (i > 0) {
                result += ", ";
            }
            result += arguments[i]->toString();
        }
        result += ")";
        return result;
    }

    NodeType getType() const override {
        return NodeType::FUNCTION;
    }

    std::shared_ptr<ExpressionNode> clone() const override {
        std::vector<std::shared_ptr<ExpressionNode>> clonedArgs;
        clonedArgs.reserve(arguments.size());
        for (const auto& arg : arguments) {
            clonedArgs.push_back(arg->clone());
        }
        return std::make_shared<FunctionNode>(name, std::move(clonedArgs));
    }

    std::shared_ptr<ExpressionNode> simplify() override {
        // Simplify all arguments
        std::vector<std::shared_ptr<ExpressionNode>> simplifiedArgs;
        simplifiedArgs.reserve(arguments.size());
        bool allArgsConstant = true;
        
        for (const auto& arg : arguments) {
            auto simplifiedArg = arg->simplify();
            simplifiedArgs.push_back(simplifiedArg);
            if (!simplifiedArg->isConstant()) {
                allArgsConstant = false;
            }
        }

        // If all arguments are constants, evaluate the function
        if (allArgsConstant) {
            try {
                auto functionNode = std::make_shared<FunctionNode>(name, simplifiedArgs);
                double result = functionNode->evaluate();
                return std::make_shared<ConstantNode>(result);
            } catch (const std::exception&) {
                // If evaluation fails, return the simplified function node
            }
        }

        // Handle special function cases
        if (name == "sin") {
            if (simplifiedArgs[0]->isConstant() && simplifiedArgs[0]->getConstantValue() == 0) {
                return std::make_shared<ConstantNode>(0); // sin(0) = 0
            }
        } else if (name == "cos") {
            if (simplifiedArgs[0]->isConstant() && simplifiedArgs[0]->getConstantValue() == 0) {
                return std::make_shared<ConstantNode>(1); // cos(0) = 1
            }
        } else if (name == "sqrt") {
            if (simplifiedArgs[0]->isConstant()) {
                double val = simplifiedArgs[0]->getConstantValue();
                if (val == 0 || val == 1) {
                    return std::make_shared<ConstantNode>(val); // sqrt(0) = 0, sqrt(1) = 1
                }
            }
        }

        return std::make_shared<FunctionNode>(name, std::move(simplifiedArgs));
    }

    const std::string& getName() const {
        return name;
    }

    const std::vector<std::shared_ptr<ExpressionNode>>& getArguments() const {
        return arguments;
    }
};

// Token types for the parser
enum class TokenType {
    NUMBER,
    VARIABLE,
    OPERATOR,
    FUNCTION,
    LEFT_PAREN,
    RIGHT_PAREN,
    COMMA,
    END
};

struct Token {
    TokenType type;
    std::string value;
};

// Expression parser class
class ExpressionParser {
private:
    std::string expression;
    size_t position;
    Token currentToken;

    // Tokenize the expression
    void getNextToken() {
        while (position < expression.size() && std::isspace(expression[position])) {
            ++position;
        }

        if (position >= expression.size()) {
            currentToken = {TokenType::END, ""};
            return;
        }

        char c = expression[position];

        if (std::isdigit(c) || c == '.') {
            std::string number;
            while (position < expression.size() && 
                  (std::isdigit(expression[position]) || expression[position] == '.')) {
                number += expression[position++];
            }
            currentToken = {TokenType::NUMBER, number};
            return;
        }

        if (std::isalpha(c)) {
            std::string identifier;
            while (position < expression.size() && 
                  (std::isalnum(expression[position]) || expression[position] == '_')) {
                identifier += expression[position++];
            }
            
            if (position < expression.size() && expression[position] == '(') {
                currentToken = {TokenType::FUNCTION, identifier};
            } else {
                currentToken = {TokenType::VARIABLE, identifier};
            }
            return;
        }

        if (c == '+' || c == '-' || c == '*' || c == '/' || c == '^') {
            currentToken = {TokenType::OPERATOR, std::string(1, c)};
            ++position;
            return;
        }

        if (c == '(') {
            currentToken = {TokenType::LEFT_PAREN, "("};
            ++position;
            return;
        }

        if (c == ')') {
            currentToken = {TokenType::RIGHT_PAREN, ")"};
            ++position;
            return;
        }

        if (c == ',') {
            currentToken = {TokenType::COMMA, ","};
            ++position;
            return;
        }

        throw std::runtime_error("Invalid character in expression: " + std::string(1, c));
    }

    // Parse an expression
    std::shared_ptr<ExpressionNode> parseExpression() {
        auto node = parseTerm();

        while (currentToken.type == TokenType::OPERATOR && 
              (currentToken.value == "+" || currentToken.value == "-")) {
            char op = currentToken.value[0];
            getNextToken();
            auto right = parseTerm();
            node = std::make_shared<BinaryOperatorNode>(op, node, right);
        }

        return node;
    }

    // Parse a term
    std::shared_ptr<ExpressionNode> parseTerm() {
        auto node = parseFactor();

        while (currentToken.type == TokenType::OPERATOR && 
              (currentToken.value == "*" || currentToken.value == "/")) {
            char op = currentToken.value[0];
            getNextToken();
            auto right = parseFactor();
            node = std::make_shared<BinaryOperatorNode>(op, node, right);
        }

        return node;
    }

    // Parse a factor
    std::shared_ptr<ExpressionNode> parseFactor() {
        auto node = parsePrimary();

        while (currentToken.type == TokenType::OPERATOR && currentToken.value == "^") {
            getNextToken();
            auto right = parsePrimary();
            node = std::make_shared<BinaryOperatorNode>('^', node, right);
        }

        return node;
    }

    // Parse a primary expression
    std::shared_ptr<ExpressionNode> parsePrimary() {
        if (currentToken.type == TokenType::NUMBER) {
            double value = std::stod(currentToken.value);
            getNextToken();
            return std::make_shared<ConstantNode>(value);
        }

        if (currentToken.type == TokenType::VARIABLE) {
            std::string name = currentToken.value;
            getNextToken();
            return std::make_shared<VariableNode>(name);
        }

        if (currentToken.type == TokenType::FUNCTION) {
            std::string name = currentToken.value;
            getNextToken();  // consume function name

            if (currentToken.type != TokenType::LEFT_PAREN) {
                throw std::runtime_error("Expected '(' after function name");
            }
            getNextToken();  // consume '('

            std::vector<std::shared_ptr<ExpressionNode>> args;
            if (currentToken.type != TokenType::RIGHT_PAREN) {
                args.push_back(parseExpression());
                while (currentToken.type == TokenType::COMMA) {
                    getNextToken();  // consume ','
                    args.push_back(parseExpression());
                }
            }

            if (currentToken.type != TokenType::RIGHT_PAREN) {
                throw std::runtime_error("Expected ')' after function arguments");
            }
            getNextToken();  // consume ')'

            return std::make_shared<FunctionNode>(name, std::move(args));
        }

        if (currentToken.type == TokenType::OPERATOR && 
           (currentToken.value == "+" || currentToken.value == "-")) {
            char op = currentToken.value[0];
            getNextToken();
            auto operand = parsePrimary();
            return std::make_shared<UnaryOperatorNode>(op, operand);
        }

        if (currentToken.type == TokenType::LEFT_PAREN) {
            getNextToken();  // consume '('
            auto node = parseExpression();
            
            if (currentToken.type != TokenType::RIGHT_PAREN) {
                throw std::runtime_error("Expected ')'");
            }
            getNextToken();  // consume ')'
            
            return node;
        }

        throw std::runtime_error("Unexpected token: " + currentToken.value);
    }

public:
    ExpressionParser(const std::string& expr) : expression(expr), position(0) {
        getNextToken();
    }

    std::shared_ptr<ExpressionNode> parse() {
        auto node = parseExpression();
        
        if (currentToken.type != TokenType::END) {
            throw std::runtime_error("Unexpected token at end of expression: " + currentToken.value);
        }
        
        return node;
    }
};

// Evaluator class
class Evaluator {
private:
    std::shared_ptr<ExpressionNode> root;

public:
    Evaluator(const std::string& expression) {
        try {
            ExpressionParser parser(expression);
            root = parser.parse();
        } catch (const std::exception& e) {
            throw std::runtime_error("Parser error: " + std::string(e.what()));
        }
    }

    double evaluate() const {
        try {
            return root->evaluate();
        } catch (const std::exception& e) {
            throw std::runtime_error("Evaluation error: " + std::string(e.what()));
        }
    }

    std::string toString() const {
        return root->toString();
    }

    void simplify() {
        root = root->simplify();
    }

    std::string visualize() const {
        std::stringstream result;
        visualizeNode(root, "", true, result);
        return result.str();
    }

private:
    void visualizeNode(const std::shared_ptr<ExpressionNode>& node, 
                       const std::string& prefix, bool isLast, 
                       std::stringstream& result) const {
        result << prefix;
        result << (isLast ? "└── " : "├── ");

        switch (node->getType()) {
            case NodeType::CONSTANT: {
                auto constantNode = std::dynamic_pointer_cast<ConstantNode>(node);
                result << "Constant: " << constantNode->toString() << "\n";
                break;
            }
            case NodeType::VARIABLE: {
                auto variableNode = std::dynamic_pointer_cast<VariableNode>(node);
                result << "Variable: " << variableNode->toString() << "\n";
                break;
            }
            case NodeType::UNARY_OPERATOR: {
                auto unaryNode = std::dynamic_pointer_cast<UnaryOperatorNode>(node);
                result << "Unary Op: " << unaryNode->getOperator() << "\n";
                visualizeNode(unaryNode->getOperand(), 
                              prefix + (isLast ? "    " : "│   "), true, result);
                break;
            }
            case NodeType::BINARY_OPERATOR: {
                auto binaryNode = std::dynamic_pointer_cast<BinaryOperatorNode>(node);
                result << "Binary Op: " << binaryNode->getOperator() << "\n";
                visualizeNode(binaryNode->getLeft(), 
                              prefix + (isLast ? "    " : "│   "), false, result);
                visualizeNode(binaryNode->getRight(), 
                              prefix + (isLast ? "    " : "│   "), true, result);
                break;
            }
            case NodeType::FUNCTION: {
                auto funcNode = std::dynamic_pointer_cast<FunctionNode>(node);
                result << "Function: " << funcNode->getName() << "\n";
                const auto& args = funcNode->getArguments();
                for (size_t i = 0; i < args.size(); ++i) {
                    visualizeNode(args[i], 
                                  prefix + (isLast ? "    " : "│   "), 
                                  i == args.size() - 1, result);
                }
                break;
            }
        }
    }

public:
    std::string traceEvaluation() const {
        std::stringstream result;
        result << "Evaluation trace for: " << root->toString() << "\n";
        result << "-------------------------------------------\n";
        traceEvaluationNode(root, 0, result);
        result << "-------------------------------------------\n";
        result << "Final result: " << root->evaluate() << "\n";
        return result.str();
    }

private:
    double traceEvaluationNode(const std::shared_ptr<ExpressionNode>& node, 
                              int depth, std::stringstream& result) const {
        std::string indent(depth * 2, ' ');
        
        switch (node->getType()) {
            case NodeType::CONSTANT: {
                auto constantNode = std::dynamic_pointer_cast<ConstantNode>(node);
                double value = constantNode->evaluate();
                result << indent << "Constant " << constantNode->toString() 
                      << " = " << value << "\n";
                return value;
            }
            case NodeType::VARIABLE: {
                auto variableNode = std::dynamic_pointer_cast<VariableNode>(node);
                double value = variableNode->evaluate();
                result << indent << "Variable " << variableNode->toString() 
                      << " = " << value << "\n";
                return value;
            }
            case NodeType::UNARY_OPERATOR: {
                auto unaryNode = std::dynamic_pointer_cast<UnaryOperatorNode>(node);
                result << indent << "Evaluating unary operator " 
                      << unaryNode->getOperator() << " on:\n";
                double operandValue = traceEvaluationNode(unaryNode->getOperand(), depth + 1, result);
                double value = unaryNode->evaluate();
                result << indent << "Result of " << unaryNode->getOperator() 
                      << "(" << operandValue << ") = " << value << "\n";
                return value;
            }
            case NodeType::BINARY_OPERATOR: {
                auto binaryNode = std::dynamic_pointer_cast<BinaryOperatorNode>(node);
                result << indent << "Evaluating binary operator " 
                      << binaryNode->getOperator() << " on:\n";
                double leftValue = traceEvaluationNode(binaryNode->getLeft(), depth + 1, result);
                double rightValue = traceEvaluationNode(binaryNode->getRight(), depth + 1, result);
                double value = binaryNode->evaluate();
                result << indent << "Result of " << leftValue << " " 
                      << binaryNode->getOperator() << " " << rightValue 
                      << " = " << value << "\n";
                return value;
            }
            case NodeType::FUNCTION: {
                auto funcNode = std::dynamic_pointer_cast<FunctionNode>(node);
                result << indent << "Evaluating function " << funcNode->getName() << " with arguments:\n";
                std::vector<double> argValues;
                for (const auto& arg : funcNode->getArguments()) {
                    argValues.push_back(traceEvaluationNode(arg, depth + 1, result));
                }
                double value = funcNode->evaluate();
                result << indent << "Result of " << funcNode->getName() << "(";
                for (size_t i = 0; i < argValues.size(); ++i) {
                    if (i > 0) result << ", ";
                    result << argValues[i];
                }
                result << ") = " << value << "\n";
                return value;
            }
        }
        return 0.0;  // Should never reach here
    }
};