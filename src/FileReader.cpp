
#include "FileReader.hpp"

#include <string>

// 初始化参数<string, int> -> <string, string>
void FileReader::registerIntParameter(const std::string &key, int init)                       // 将 int 转化成 String
{
    std::stringstream ss;                                                                     // 类型转换的中介
    ss << init;                                                                               // 将整数 init 写入 ss 中
    parameters[key] = ss.str();
}


void FileReader::registerdoubleParameter(const std::string &key, double init)
{
    std::stringstream ss;
    ss << init;
    parameters[key] = ss.str();
}


void FileReader::registerStringParameter(const std::string &key, const std::string &init)     
{
    parameters[key] = init;
}

// 设置参数
void FileReader::setParameter(const std::string &key, const std::string &in)
{

    parameters[key] = in;
}


void FileReader::setParameter(const std::string &key, double in)
{
    std::stringstream ss;
    ss << in;

    parameters[key] = ss.str();
}


void FileReader::setParameter(const std::string &key, int in)
{
    std::stringstream ss;
    ss << in;
    parameters[key] = ss.str();
}

// 读取指定配置文件，解析其中的键值对将其加载到程序中
bool FileReader::readFile(const std::string &name)
{
    
    std::fstream fileInput(name, std::ios_base::in);                                     // 创建一个文件流对象，以读取模式 in 打开 name 文件
    std::string readout, elem1, elem2;                                                   // readout 用来暂存从文件中读取的一行字符串，elem1放键，elem2放值

    while (fileInput.is_open())                                                          // 文件为打开状态
    {
        while (std::getline(fileInput, readout))                                         // 还能读取到一行字符串
        {
            std::istringstream ss(readout, std::istringstream::in);                      // 将 readout 变成流
            ss >> elem1 >> elem2;                                                        // 将第一个值给 elem1，第二个值给 elem2
            if(!elem1.empty() && !elem2.empty() && elem1[0] != '#' and elem2[0] != '#')  // 注释及空行检查
            {
                try
                {
                    double i = std::stod(elem2);                                         // string -> double
                    registerdoubleParameter(elem1, i);                                   // 赋值
                    setParameter(elem1, i);
                }
                catch(...)
                {
                    registerStringParameter(elem1, elem2);
                    setParameter(elem1,elem2);
                }
            }
        }
        fileInput.close();
    }

    return true;
}

// 遍历打印出当前 FileReader 类中存储的所有键值对，常用于调试
void FileReader::printParameters() const
{
    for (auto it = parameters.begin(); it != parameters.end(); it++)                           // auto 让编译器自动推断类型
    {
        std::cout << it->first << '\t' << it->second << std::endl;
    }
}

bool FileReader::findparameters(std::string param)                                            // 检查指定的 key 是否存在于当前 map 中
{

    auto it = parameters.find(param);                                                         // 找到了指向该元素，没找到指向 parameters.end()
    if(parameters.find(param) == parameters.end())
        return false;
    else
        return true;
}

std::string FileReader::findstring(std::string param) const
{
    if(parameters.find(param) == parameters.end())
        return "noslip";
    else
        return getStringParameter(param);
}