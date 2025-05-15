#include "express.cpp"  
#include <iostream>
#include <string>
#include <map>
#include <memory>
#include <limits>

// Rest of main.cpp unchanged...
// Forward declarations are already in the original file
// Implementation of all the classes (ExpressionNode, ConstantNode, etc.) are in the original file

// Simple command-line interface for the expression evaluator
int main() {
    // Initialize built-in functions
    initFunctions();
    
    std::cout << "Expression Tree Evaluator" << std::endl;
    std::cout << "=========================" << std::endl;
    std::cout << "Type 'help' for command list or 'exit' to quit" << std::endl;
    
    std::string input;
    
    while (true) {
        std::cout << "\n> ";
        std::getline(std::cin, input);
        
        if (input == "exit") {
            break;
        } else if (input == "help") {
            std::cout << "Available commands:" << std::endl;
            std::cout << "  evaluate <expression>  - Evaluate a mathematical expression" << std::endl;
            std::cout << "  set <var> <value>      - Set a variable value" << std::endl;
            std::cout << "  vars                   - List all defined variables" << std::endl;
            std::cout << "  functions              - List all available functions" << std::endl;
            std::cout << "  simplify <expression>  - Simplify an expression" << std::endl;
            std::cout << "  visualize <expression> - Visualize expression tree" << std::endl;
            std::cout << "  trace <expression>     - Show evaluation trace" << std::endl;
            std::cout << "  exit                   - Exit the program" << std::endl;
        } else if (input.substr(0, 8) == "evaluate") {
            if (input.length() <= 9) {
                std::cout << "Error: Missing expression" << std::endl;
                continue;
            }
            
            std::string expr = input.substr(9);
            try {
                Evaluator evaluator(expr);
                double result = evaluator.evaluate();
                std::cout << "Result: " << result << std::endl;
            } catch (const std::exception& e) {
                std::cout << "Error: " << e.what() << std::endl;
            }
        } else if (input.substr(0, 3) == "set") {
            size_t pos1 = input.find(' ', 4);
            if (pos1 == std::string::npos) {
                std::cout << "Error: Invalid format. Use 'set <var> <value>'" << std::endl;
                continue;
            }
            
            std::string var = input.substr(4, pos1 - 4);
            std::string valueStr = input.substr(pos1 + 1);
            
            try {
                double value = std::stod(valueStr);
                variables[var] = value;
                std::cout << "Variable '" << var << "' set to " << value << std::endl;
            } catch (const std::exception& e) {
                std::cout << "Error: Invalid value" << std::endl;
            }
        } else if (input == "vars") {
            if (variables.empty()) {
                std::cout << "No variables defined" << std::endl;
            } else {
                std::cout << "Defined variables:" << std::endl;
                for (const auto& [name, value] : variables) {
                    std::cout << "  " << name << " = " << value << std::endl;
                }
            }
        } else if (input == "functions") {
            std::cout << "Available functions:" << std::endl;
            for (const auto& [name, funcInfo] : functions) {
                std::cout << "  " << name << " (" << funcInfo.second << " args)" << std::endl;
            }
        } else if (input.substr(0, 8) == "simplify") {
            if (input.length() <= 9) {
                std::cout << "Error: Missing expression" << std::endl;
                continue;
            }
            
            std::string expr = input.substr(9);
            try {
                Evaluator evaluator(expr);
                std::cout << "Original: " << evaluator.toString() << std::endl;
                evaluator.simplify();
                std::cout << "Simplified: " << evaluator.toString() << std::endl;
            } catch (const std::exception& e) {
                std::cout << "Error: " << e.what() << std::endl;
            }
        } else if (input.substr(0, 9) == "visualize") {
            if (input.length() <= 10) {
                std::cout << "Error: Missing expression" << std::endl;
                continue;
            }
            
            std::string expr = input.substr(10);
            try {
                Evaluator evaluator(expr);
                std::cout << evaluator.visualize() << std::endl;
            } catch (const std::exception& e) {
                std::cout << "Error: " << e.what() << std::endl;
            }
        } else if (input.substr(0, 5) == "trace") {
            if (input.length() <= 6) {
                std::cout << "Error: Missing expression" << std::endl;
                continue;
            }
            
            std::string expr = input.substr(6);
            try {
                Evaluator evaluator(expr);
                std::cout << evaluator.traceEvaluation() << std::endl;
            } catch (const std::exception& e) {
                std::cout << "Error: " << e.what() << std::endl;
            }
        } else {
            std::cout << "Unknown command. Type 'help' for available commands." << std::endl;
        }
    }
    
    return 0;
}