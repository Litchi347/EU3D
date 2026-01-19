#pragma once
#include <assert.h>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <utility>
#include <tuple>
namespace ARRAY                                                                // 将多维索引映射到一维连续内存空间中
{
    template <class Lam, size_t... Is>                                         // 模版参数，一个名为 Lam 的类型和无符号整型的多个数字，按照数字个数展开代码
    inline constexpr void Loop(Lam lam,std::index_sequence<Is...>)             // lam 相当于循环体，Is 相当于循环的索引，在编译期间循环
    {
        (lam(Is), ...);
    };

    template <class T, size_t Dim>                                             // 模版声明
    class Array
    {
    private:
        T *buf = nullptr;                                                      // 指向连续内存块的指针
        size_t Len = 0;                                                        // 数组总元素的个数
        size_t EachDim[Dim] = {};                                              // 每一维具体长度

        template <size_t... Is>
        inline const void PrintEachDim(std::index_sequence<Is...>)
        {
            auto PrintDim = [](auto i)
            { std::cout << i << '\t'; };
            (PrintDim(EachDim[Is]), ...);
        }
        // 行优先映射一维地址
        inline const size_t CalIndex(std::array<int, Dim> &CartIndex)
        {
            size_t FlatIndex = 0;
            Loop(
                [&](auto i) constexpr __attribute__((always_inline)) {         // [&] 以引用方式捕获外部变量以便在内部修改它们，auto i 是循环当前索引，constexpr 尽可能在编译阶段完成
                    FlatIndex = (CartIndex[i] + FlatIndex) * EachDim[i + 1];
                },
                std::make_index_sequence<Dim - 1>{});                          // 生成从 0 到 Dim-2 的整数序列
            FlatIndex += CartIndex[Dim - 1];
            return FlatIndex;
        }

    public:
        Array() = default;                                                     // 默认构造函数


        template <class... EACHDIM>
        Array(EACHDIM... each)                                                 // 变长参数构造函数，根据传入的任意数量的维度大小分配一块连续内存并记录每一维的长度
        {
            constexpr bool IsInteger = 
                (true && ... && std::is_integral_v<decltype(each)>);           // 检查 each 中每一个参数是不是都是 int 类型

            static_assert(IsInteger == true);                                  // 确保 each 中每一个参数都是 int 类型
            
            static_assert(sizeof...(each) == Dim);                             // 确保传入参数类型等于 Dim


            auto GetBufLen = [&]() -> size_t
            { return(1 * ... * each); };                                       // 将所有维度长度乘起来
            Len = GetBufLen();
            // 申请连续内存空间
            buf = new T[Len];


            std::array<int, sizeof...(each)> EACH = {each...};                 // 初始化成一个 array 便于后续使用索引访问
            Loop(
                [&](auto i) constexpr
                __attribute__((always_inline)){ EachDim[i] = EACH[i]; },
                std::make_index_sequence<sizeof...(each)>{});
        }

        // 初始化
        template <class... EACHDIM>
        void Initial(EACHDIM... each)
        {
            if(buf)
            {
                delete[] buf;
                buf = nullptr;
            }

            constexpr bool IsInteger = 
                (true && ... && std::is_integral_v<decltype(each)>);

            static_assert(IsInteger == true);

            static_assert(sizeof...(each) == Dim);


            auto GetBufLen = [&]() -> size_t
            { return (1 * ... * each); };
            Len = GetBufLen();

            buf = new T[Len];


            std::array<int, sizeof...(each)> EACH = {each...};
            Loop(
                [&](auto i) constexpr
                __attribute__((always_inline)) { EachDim[i] = EACH[i]; },
                std::make_index_sequence<sizeof...(each)>{});
        }


        template <class T_t, size_t Dim_t>
        Array(const Array<T_t, Dim_t> &array)                                  // 允许不同数据类型但维度相同的数组进行转化拷贝
        {

            constexpr auto IsConvert = std::is_convertible_v<T_t, T>;
            static_assert(IsConvert);                                          // 确保源类型可以安全转换成目标类型

            static_assert(Dim == Dim_t);

            Len = array.GetLen();
            buf = new T[Len];
            auto array_ptr = array.Getbuf();                                   // 源数组底层的原始指针
            for (int i = 0; i < Len; i++)
                buf[i] = array_ptr[i];                                         // 自动发生隐式类型转换
        }








        // 析构函数，在对象生命周期结束时被自动调用
        ~Array()
        {
            if(buf)
                delete[] buf;
            // std::cout << "Objects is deleted\n";
        }

        void Delete()
        {
            if(buf)
                delete[] buf;
            // std::cout << "Deleted\n";
        }

        // 打印维度信息
        void GetDim()
        {
            std::cout << "Dim:" << Dim << "\nEach Dim:";
            PrintEachDim(std::make_index_sequence<Dim>{});
            std::cout << std::endl;
        }


        size_t GetLen() const{ return Len; }


        T *Getbuf() const { return buf; }

        // 重载 () 运算符
        template <class... INDEX>
        T &operator()(INDEX... index)
        {
            constexpr bool IsInteger = 
                (true && ... && std::is_integral_v<decltype(index)>);

            static_assert(IsInteger == true);                                  // 确保传入的坐标是整数
            
            static_assert(sizeof...(index) == Dim);

            std::array<int, Dim> CartDim = {index...};
            auto Flat = CalIndex(CartDim);                                     // 一维索引的映射

// Debug模式下的安全防护
// 只有在 #define Debug 时才会生效
#ifdef Debug
            auto IsInRange = true;
            Loop(
                [&](auto i)
                {
                    IsInRange =                                                // 遍历输入的每一个维度坐标检查是否在 0 到 EachDim[i] 之间
                        IsInRange && CartDim[i] < EachDim[i] && CartDim[i] >= 0;
                },
                std::make_index_sequence<Dim>{});
            assert(IsInRange);
            assert(Flat >= 0 && Flat < Len);
#endif
            return buf[Flat];
        }

        // 重载 [] 运算符
        T &operator[](size_t Flat)
        {
#ifdef Debug
            assert(Flat >= 0 && Flat < Len);
#endif
            return buf[Flat];
        } 

        // 重载 = 运算符，类型转换赋值，不同数据类型
        template <class T_t>
        Array<T, Dim> &operator=(const Array<T_t, Dim> &array)
        {


            constexpr auto IsConvert = std::is_convertible_v<T_t, T>;
            static_assert(IsConvert);

            if (buf)
                delete[] buf;
            Len = array.GetLen();
            buf = new T[Len];
            for (int i = 0; i < Dim; i++)
                EachDim[i] = array.EachDim[i];
            auto array_ptr = array.Getbuf();
            for(int i = 0; i < Len; i++)
                buf[i] = array_ptr[i];
            return *this;
        }

        // 重载 = 运算符，标准拷贝赋值，相同数据类型
        Array<T, Dim> &operator=(const Array<T, Dim> &array)
        {
            // std::cout << "COPY\n";
            if (this == &array)
                return *this;
            if (buf)
                delete[] buf;
            Len = array.GetLen();
            buf = new T[Len];
            for (int i = 0 ; i < Dim; i++)
                EachDim[i] = array.EachDim[i];
            auto array_ptr = array.Getbuf();
            for (int i = 0; i < Len; i++)
                buf[i] = array_ptr[i];
            return *this;
        } 

        // 重载 == 运算符，通过逐个对比内存中的元素，判断两个数组的内容是否完全一致
        template <class T_t,size_t Dim_t>
        bool operator==(Array<T_t, Dim_t> &array)
        {

            constexpr auto IsSame = std::is_same_v<T, T_t>;
            static_assert(IsSame);
            // 确保维度相同
            static_assert(Dim == Dim_t);

            if(Len != array.GetLen())
                return false;
            else
            {
                auto result = true;
                auto array_ptr = array.Getbuf();
                for (int i = 0; i < Len; i++)
                    result = result && (buf[i] == array_ptr[i]);
                return result;
            }
        }

        // 使用一个特定的数值填充 array
        template <class Type>
        void Fill(Type value)
        {
            constexpr auto TypeOK = std::is_convertible_v<Type, T>;
            static_assert(TypeOK);
            for (int i = 0; i < Len; i++)
                buf[i] = value;
        }

        // 找最大值
        inline const T MaxValue() const
        {
            T maxValue = buf[0];
            for (size_t i = 1; i < Len; i++)
            {
                if(buf[i] > maxValue)
                {
                    maxValue = buf[i];
                }
            }
            return maxValue;
        }
        // 找最大值所在的位置
        inline const size_t MaxPosition() const
        {
            T maxValue = buf[0];
            size_t maxPosition = 0;
            for (size_t i = 1; i < Len; i ++)
            {
                if(buf[i] > maxValue)
                {
                    maxValue = buf[i];
                    maxPosition = i;
                }
            }
            return maxPosition;
        }


        inline const T MinValue() const
        {
            T minValue = buf[0];
            for (size_t i = 1; i < Len; i++)
            {
                if(buf[i] < minValue)
                {
                    minValue = buf[i];
                }
            }
            return minValue;
        }


        inline const T AveValue() const
        {
            return Sum() / Len;

        }

        inline const T Sum() const
        {
            T sum = 0;
            for (size_t i = 0; i < Len; i++)
            {
                sum += buf[i];
            }
            return sum;
        }

        inline const T SumNoBoundary(int bc) const                             // 去除边缘层后内部数据部分化学计算步数总和
        {
            T sum = 0;
            for (size_t i = 0; i < Len; i++)
            {
                size_t z = i % EachDim[2];                                     // z 坐标
                size_t y = (i / EachDim[2]) % EachDim[1];                      // y 坐标
                size_t x = i / (EachDim[2] * EachDim[1]);                      // x 坐标
                if (x >= bc && y >= bc && z >= bc && 
                    x < EachDim[0] - bc &&
                    y < EachDim[1] - bc &&
                    z < EachDim[2] - bc)
                    sum += buf[i];
            }
            return sum;
        }


        inline const T SumPositive() const
        {
            T sump = 0;
            for (size_t i = 0; i < Len; i++)
            {
                if (buf[i] > 0)
                {
                    sump += buf[i];
                }
            }
            return sump;
        }

        // 判断有无空值
        inline const bool IsNan() const
        {
            for (size_t i = 0; i < Len; i++)
            {
                if(std::isnan(buf[i]))
                {
#ifdef Debug
                    std::cout << "Exisr NAN\n";
#endif
                    return true;
                }
            }
#ifdef Debug
            std::cout << "No NAN\n";
#endif
            return false;
        }


        inline void Print()
        {
            for (size_t i = 0; i < Len; i++)
                std::cout << buf[i] << '\t';
            std::cout << '\n';
        }
        inline size_t GetSize() const { return Len; }
        inline size_t GetDim() const { return Dim; }
        inline size_t GetNi() const { return EachDim[0]; }
        inline size_t GetNj() const { return EachDim[1]; }
        inline size_t GetNk() const { return EachDim[2]; }
        std::tuple<size_t, size_t, size_t> Get3DIndices(size_t m)              // 得到三维坐标
        {
            size_t k = m % EachDim[2];
            size_t j = (m / EachDim[2]) % EachDim[1];
            size_t i = m / (EachDim[2] * EachDim[1]);
            return std::make_tuple(i, j, k);
        }
    };
}