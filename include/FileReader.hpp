
















#ifndef FILEREADER_HH                                                                    // 头文件保护符
#define FILEREADER_HH                                                                    // 防止该头文件在一个编译单元中被包含多次，避免重复定义错误





# include <string>
# include <cstring>
# include <iostream>
# include <fstream>
# include <map>
# include <sstream>

class FileReader                                                                         // 定义了一个FileReader 类用来读取输入文件
{
public:                                                                                  // 可以被外部访问
    // 初始化参数
    void registerIntParameter(const std::string &key, int init);


    void registerdoubleParameter(const std::string &key, double init);


    void registerStringParameter(const std::string &key, const std::string &init = "");

    // 设置参数
    void setParameter(const std::string &key, const std::string &in);


    void setParameter(const std::string &key,double in);


    void setParameter(const std::string &key,int in);

    // 获取参数
    inline int getIntParameter(const std::string &key) const;                            // inline 表示在头文件中直接定义


    inline double getdoubleParameter(const std::string &key) const;                      // const 表示在函数执行过程中不会改变类的成员变量


    inline std::string getStringParameter(const std::string &key) const;


    bool readFile(const std::string &name);


    void printParameter() const;
    bool findparameters(std::string param);
    std::string findstring(std::string param) const;

private:
    std::map<std::string, std::string> parameters;
};

inline int FileReader::getIntParameter(const std::string &key) const
{
    auto it = parameters.find(key);                                                      // 在 parameters 的 map 中搜索键为 key 的项 it->first 获取键， it->second 获取值 
    try
    {
        std::stoi(it -> second);                                                         // 将其对应的值转换成 int 类型
        std::cout << it -> first << '\t' << it -> second << std::endl;
    }
    catch(...)
    {
        return 0;
    }

    return std::stoi(it -> second);
}

inline double FileReader::getdoubleParameter(const std::string &key) const
{
    auto it = parameters.find(key);
    try
    {
        std::stoi(it -> second);
        std::cout << it -> first << '\t' << it -> second << std::endl;
    }
    catch(...)
    {
        return 0;
    }

    return std::stod(it -> second);
}

inline std::string FileReader::getStringParameter(const std::string &key) const
{
    auto it = parameters.find(key);
    std::cout << it -> first << '\t' << it -> second << std::endl;
    return it -> second;
}

#endif
