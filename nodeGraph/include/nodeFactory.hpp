#pragma once
#include "NodeGraph.hpp"
#include <functional>
// 工厂类 - 负责类的创建
class NodeFactory {
private:
    // 定义创建函数的类型
    using CreateFunction = std::function<Node* ()>;
    // 存储类名和创建函数的映射
    static std::unordered_map<std::string, CreateFunction> creators;

public:
    // 注册类的创建函数
    static void registerClass(const std::string& className, CreateFunction func);
    // 根据类名创建对象
    static Node* createObject(const std::string& className);
};