#include "nodeFactory.hpp"
#include "InterpolateCurve.hpp"
#include "MergeNode.hpp"
#include "numberNode.hpp"
#include "pointNode.hpp"


// 注册类的创建函数
void NodeFactory::registerClass(const std::string& className, CreateFunction func) 
{
	creators[className] = func;
}

// 根据类名创建对象
Node* NodeFactory::createObject(const std::string& className)
{
	auto it = creators.find(className);
	if (it != creators.end()) {
		return it->second();
	}
	std::cerr << "Error: Class " << className << " not registered!" << std::endl;
	return nullptr;
}



// 初始化静态成员
std::unordered_map<std::string, NodeFactory::CreateFunction> NodeFactory::creators;

// 注册函数模板
template<typename T>
Node* createInstance() 
{
    return new T(std::string(""));
}

// 注册宏 - 简化注册过程
#define REGISTER_CLASS(className) \
    static struct className##Registrar { \
        className##Registrar() { \
            NodeFactory::registerClass(#className, &createInstance<className>); \
        } \
    } className##RegistrarInstance;

// 注册具体的类
REGISTER_CLASS(Ellipse)
REGISTER_CLASS(MergeNode)
REGISTER_CLASS(ListItem)
REGISTER_CLASS(Rotate)
REGISTER_CLASS(Scale)
REGISTER_CLASS(CurveSplit)
REGISTER_CLASS(Loft)
REGISTER_CLASS(Sweep)
REGISTER_CLASS(Flatten)
REGISTER_CLASS(Line)
REGISTER_CLASS(InterpolateCurve)
REGISTER_CLASS(CurveIntCurve)
REGISTER_CLASS(NumberNode)
REGISTER_CLASS(Ranges)
REGISTER_CLASS(Extrude)
REGISTER_CLASS(ConstructPoint)
REGISTER_CLASS(DistancePt)
REGISTER_CLASS(EndsPoints)
REGISTER_CLASS(PrepFrame)
REGISTER_CLASS(Rectangle)
