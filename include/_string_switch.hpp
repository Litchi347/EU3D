#include <iostream>

using namespace std;

typedef std::uint64_t hash_t;                                                            // 将 64 为无符号整数定义为哈希类型

constexpr hash_t prime = 0x100000001B3ull;                                               // 素数
constexpr hash_t basis = 0xCBF29CE484222325ull;                                          // 初始偏移量
// 运行期哈希函数
hash_t hash_(char const *str)                                                            // 计算字符的哈希值
{
    hash_t ret{basis};                                                                   // 初始值为基准值

    while (*str)                                                                          // 遍历字符串直至遇到结束符'\0'
    {
        ret ^= *str;                                                                     // 将当前字符与结果进行异或
        ret *= prime;                                                                    // 将结果乘以素数，进行“混淆”
        str ++;                                                                          // 指针后移
    }

    return ret;
}
// 编译期哈希函数，字符串直接变为数字，运行时 0 开销
constexpr hash_t hash_compile_time(char const *str, hash_t last_value = basis)
{
    return *str ? hash_compile_time(str + 1, (*str ^ last_value) * prime) : last_value;
}
// 自定义字面量运算符
constexpr unsigned long long operator"" _hash(char const *p, size_t)
{
    return hash_compile_time(p);
}
