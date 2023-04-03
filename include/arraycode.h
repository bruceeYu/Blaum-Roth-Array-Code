
#ifndef __ARRAYCODE_H
#define __ARRAYCODE_H

#include</home/brucee/Public/ISA-L/include/raid.h>
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdbool.h>

/// @brief 计算数学意义上的余数,例如-5模7等于+2
/// @param value 整数value
/// @param prime 整数prime
/// @return value对prime的余数
int64_t mod(int64_t value, int64_t prime);

/// @brief 取大值 
/// @param a 整数a
/// @param b 整数b
/// @return 整数a和整数b的较大值
int64_t max(int64_t a, int64_t b);

/// @brief 取小值
/// @param a 整数a
/// @param b 整数b
/// @return 整数a和整数b的较小值
int64_t min(int64_t a, int64_t b);

/// @brief 整数a和整数b距离
/// @param a 整数a
/// @param b 整数b
/// @return 整数a和整数b差的绝对值
int64_t diff(int64_t a, int64_t b);

/// @brief 交换整数a和整数b的值
/// @param pa 指向整数a的指针
/// @param pb 指向整数b的指针
void swap(int64_t* pa, int64_t* pb);

/// @brief 计算多项式ptr第index个符号的首地址
/// @param ptr 多项式，步长为1
/// @param index 第index个符号
/// @param prime index需对prime取模
/// @param packet_size 符号大小
/// @return 第index个符号的首地址
uint8_t* locate_packet(const uint8_t* ptr, int64_t index, int64_t prime, int64_t packet_size);

/// @brief 对两个非零项的因式乘积的每个因式提取公共项，最终化简为x^pairs[0](1+x^pairs[1])...(1+x^pairs[len-1])
/// @param pairs 对应原始因式乘积(x^pairs[0]+x^pairs[1])(x^pairs[0]+x^pairs[1])...(x^pairs[0]+x^pairs[len-1])
/// @param len 指数数组pairs的长度或元素个数
/// @param prime 用于对指数取模，例如(x^3+x^2)(x^3+x^5)(x^3+x^6)模(1+x^7)提取公共项后为x(1+x)(1+x^2)(1+x^3)
void extract_commom_factor(int64_t* pairs, const int64_t len, const int64_t prime);

/// @brief 对两个非零项的因式乘积进行合并，最终化简为x^pairs[0](1+x^pairs[1])...
/// @param pairs 对应原始因式乘积(x^pairs[0]+x^pairs[1])(x^pairs[0]+x^pairs[1])...(x^pairs[0]+x^pairs[len-1])
/// @param len 指数数组pairs的长度或元素个数
/// @param prime 用于对指数取模
/// @return 最简因式乘积所对应的指数的个数，例如x^4(1+x^3)，指数为4和3,有两个指数
int64_t integrate(int64_t* pairs, const int64_t len, const int64_t prime);

/// @brief 实例化一个含有n个码元的码字，其中冗余码元个数为redundancy，每个码元包含prime个大小为packet_size的符号
/// @param prime 一个码元包含的符号个数
/// @param n 码元个数
/// @param redundancy 冗余码元个数
/// @param packet_size 符号大小
/// @return 存放该码字的二维数组（指针）
uint8_t** create_codeword(int64_t prime, int64_t n, int64_t redundancy, int64_t packet_size);

/// @brief 释放一个含有n个码元的码字
/// @param codeword 码字
/// @param n 码元个数
/// @return 空指针
uint8_t** destroy_codeword(uint8_t** codeword, int64_t n);

/// @brief 申请amount个chunk_size字节的元素
/// @param amount 元素个数
/// @param chunk_size 元素大小
/// @return 存放amount个元素的空间指针
uint8_t** create_elements(int64_t amount, int64_t chunk_size);

/// @brief 回收amount个元素
/// @param mem 空间指针（二级指针）
/// @param amount 元素个数
/// @return 空指针
uint8_t** destroy_elements(uint8_t** mem, int64_t amount);

/// @brief 将多项式ptr指向的连续chunk_size字节的空间清零
/// @param ptr 多项式ptr
/// @param chunk_size 空间字节数
void init_chunk(uint8_t* ptr, int64_t chunk_size);

/// @brief 随机化多项式ptr前prime-1个符号，最后一个符号清零
/// @param ptr 多项式ptr
/// @param prime 符号个数
/// @param packet_size 符号大小
void rand_chunk(uint8_t* ptr, int64_t prime, int64_t packet_size);

/// @brief 将多项式ptr最后一个符号分别与前prime-1个符号进行异或
/// @param ptr 多项式ptr
/// @param prime 符号个数
/// @param packet_size 符号大小
void rectify(uint8_t* ptr, int64_t prime, int64_t packet_size);

/// @brief 将多项式ptr的全部符号循环右移shift步
/// @param ptr 多项式ptr
/// @param shift 符号移动的步长
/// @param prime 符号个数
/// @param packet_size 符号大小
void local_cyclic_shift(uint8_t* ptr, int64_t shift, int64_t prime, int64_t packet_size);

/// @brief 加等于运算 dest = dest + x^shift src 模 (1 + x^prime) 源多项式和结果多项式不能重合！！！
/// @param dest 结果多项式
/// @param src 源多项式
/// @param shift 源多项式对应符号移动的步长
/// @param prime 符号个数 
/// @param packet_size 符号大小
void add_equal(uint8_t* dest, uint8_t* src, int64_t shift, int64_t prime, int64_t packet_size);

/// @brief 除等于运算 dest = dest / (x^exp_first + x^exp_second)
/// @param dest 既是源多项式也是结果多项式
/// @param exp_first 对应 x^exp_first
/// @param exp_second 对应 x^exp_second
/// @param prime 符号个数
/// @param packet_size 符号大小
/// @param fast 值为true时，转调函数rough_divide；否则转调函数exact_divide
void divide_equal(uint8_t* dest, int64_t exp_first, int64_t exp_second, int64_t prime, int64_t packet_size, bool fast);

/// @brief 多项式dest的第i个符号为多项式src1的第i-shift1个符号和多项式src2的第i-shift2个符号的异或和
/// @param dest 结果多项式
/// @param src1 1号源多项式
/// @param shift1 1号源多项式对应符号移动的步长
/// @param src2 2号源多项式
/// @param shift2 2号源多项式对应符号移动的步长
/// @param prime 符号个数
/// @param packet_size 符号大小
void vector_add(uint8_t* dest, const uint8_t* src1, int64_t shift1, const uint8_t* src2, int64_t shift2, int64_t prime, int64_t packet_size);

/// @brief 在模1+x^prime时将乘积 fx{x^pairs[0](1+x^pairs[1])(1+x^pairs[2])...(1+x^pairs[len-1])}的结果保存于多项式gx中
/// @param gx 结果多项式
/// @param fx 源多项式
/// @param prime 符号个数
/// @param packet_size 符号大小
/// @param pairs 对应因式乘积x^pairs[0](1+x^pairs[1])(1+x^pairs[2])...(1+x^pairs[len-1])
/// @param len 因式乘积所对应的指数的个数，例如x^4(1+x^3)，指数为4和3,有两个指数
void indirect_series_mul(uint8_t* gx, const uint8_t* fx, int64_t prime, int64_t packet_size, int64_t* pairs, int64_t len);

/// @brief 在模1+x^prime时将商fx/(x^exp_first + x^exp_second)保存于多项式gx中，不保证gx含有偶数个非零项！
/// @param gx 结果多项式
/// @param fx 源多项式，要求含有偶数个非零项，不满足时使用该函数会出错！！！
/// @param prime 符号个数
/// @param packet_size 符号大小
/// @param exp_first 指数1
/// @param exp_second 指数2
void rough_divide(uint8_t* gx, const uint8_t* fx, int64_t prime, int64_t packet_size, int64_t exp_first, int64_t exp_second);

/// @brief 在模1+x^prime时将商fx/(x^exp_first + x^exp_second)保存于多项式gx中，保证gx含有偶数个非零项！
/// @param gx 结果多项式
/// @param fx 源多项式，要求含有偶数个非零项，不满足时使用该函数会出错！！！
/// @param prime 符号个数
/// @param packet_size 符号大小
/// @param exp_first 指数1
/// @param exp_second 指数2
void exact_divide(uint8_t* gx, const uint8_t* fx, int64_t prime, int64_t packet_size, int64_t exp_first, int64_t exp_second);

/// @brief 在模1+x^prime时将fx/(x^exp_first + x^exp_second)保存于多项式gx中
/// @param gx 结果多项式
/// @param fx 源多项式，要求含有偶数个非零项，不满足时使用该函数会出错！！！
/// @param prime 符号个数
/// @param packet_size 符号大小
/// @param exp_first 指数1
/// @param exp_second 指数2
/// @param fast 值为true时，转调函数rough_divide；否则转调函数exact_divide
void indirect_divide(uint8_t* gx, const uint8_t* fx, int64_t prime, int64_t packet_size, int64_t exp_first, int64_t exp_second, bool fast);

/// @brief 在模1+x^prime时将fx/(x^pairs[0](1+x^pairs[1])(1+x^pairs[2])...(1+x^pairs[len-1]))的结果就地保存在多项式fx中
/// @param fx 既是源多项式也是结果多项式，要求含有偶数个非零项，不满足时使用该函数会出错！！！函数完成时，fx=fx mod Mpx=1+x+...+x^{prime-1}
/// @param prime 符号个数
/// @param packet_size 符号大小
/// @param pairs 对应因式乘积x^pairs[0](1+x^pairs[1])(1+x^pairs[2])...(1+x^pairs[len-1])
/// @param len 因式乘积所对应的指数的个数，例如x^4(1+x^3)，指数为4和3,有两个指数
void indirect_series_divide(uint8_t* fx, int64_t prime, int64_t packet_size, int64_t* pairs, int64_t len);

/// @brief 在模Mpx=1+x+...+x^(prime-1)时将fx/(x^exp_first + x^exp_second)保存于多项式gx中
/// @param gx 结果多项式，不需要解第prime-1个符号，直接用零填充
/// @param fx 源多项式，要求第prime-1个符号等于零，不满足时使用该函数会出错！！！
/// @param prime 符号个数
/// @param packet_size 符号大小
/// @param exp_first 指数1
/// @param exp_second 指数2
void direct_divide(uint8_t* gx, const uint8_t* fx, int64_t prime, int64_t packet_size, int64_t exp_first, int64_t exp_second);

/// @brief 在模Mpx=1+x+...+x^(prime-1)时将fx/{(x^pairs[0]+x^pairs[1])(x^pairs[0]+x^pairs[1])...(x^pairs[0]+x^pairs[len-1])}的结果就地保存在多项式fx中
/// @param fx 源多项式，要求第prime-1个符号等于零，不满足时使用该函数会出错！！！
/// @param prime 符号个数
/// @param packet_size 符号大小
/// @param pairs 对应因式乘积(x^pairs[0]+x^pairs[1])(x^pairs[0]+x^pairs[1])...(x^pairs[0]+x^pairs[len-1])
/// @param len 因式乘积所对应的指数的个数，例如(x^2+x)(x^2+x^6)，指数为2，1和6,有三个指数
void direct_series_divide(uint8_t* fx, int64_t prime, int64_t packet_size, const int64_t* pairs, int64_t len);

/// @brief 在范围[0,n-1]上，把不包含在数组erasureIndex中的剩余值存入数组survivalIndex中
/// @param erasureIndex 已经选走的
/// @param lambda 数组erasureIndex的元素个数
/// @param survivalIndex 剩下的 
/// @param delta 数组survivalIndex的元素个数
void fill_survival_index(const int64_t* erasureIndex, int64_t lambda, int64_t* survivalIndex, int64_t delta);

/// @brief 计算组合数C(n,m)
/// @param n 总数 n≥m
/// @param m 取出的
/// @return 在n个数中取出m个数的组合总数
int64_t cal_combination(int64_t n, int64_t m);

/// @brief 计算码字codeword在失效列为erasure_index[0..len-1]时的amount个典型值，将典型值返回
/// @param codeword 阵列码的一个码字
/// @param prime 符号个数
/// @param n 码元个数
/// @param amount 典型值个数
/// @param packet_size 符号大小
/// @param erasure_index 失效列数组
/// @param len 失效列数组erasure_index的元素个数
/// @return 利用二维数组（指针）带回amount个典型值
uint8_t** cal_syndrome(uint8_t** codeword, int64_t prime, int64_t n, int64_t amount, int64_t packet_size, int64_t* erasure_index, int64_t len);

/// @brief 在[begin..end]中取出剩下的remain个数，将已经取出的数作为前缀，拼成一个完整的组合后存入数组patterns中
/// @param begin 范围开始
/// @param end 范围结束
/// @param amount 额定取出的个数
/// @param remain 剩余待取的个数
/// @param output 一次组合
/// @param patterns 保存所有组合的数组
/// @param seq 在结果集中属于第几个组合，编号从0开始
void simulate_aux(int64_t begin, int64_t end, int64_t amount, int64_t remain, int64_t* output, int64_t** patterns, int64_t* seq);

/// @brief 保存在[begin..end]取出amount个数的所有组合于pattern中
/// @param begin 范围开始
/// @param end 范围结束
/// @param amount 取出的
/// @param patterns 保存组合的二维数组
/// @return 组合数
int64_t simulate(int64_t begin, int64_t end, int64_t amount, int64_t** patterns);

/// @brief 判断ptr对应的内存空间是否为零
/// @param ptr 单位指针
/// @param mem_size 字节数
/// @return 内存空间全是零，返回true；否则返回false
bool equal_zero(uint8_t* ptr, int64_t mem_size);

/// @brief 计算码字的校验和
/// @param codeword 阵列码的一个码字
/// @param prime 一个码元包含的符号个数
/// @param n 码元个数
/// @param redundancy 冗余码元个数，也即典型值个数
/// @param packet_size 符号大小
/// @return 码字校验和为零，返回true；否则返回false
bool checksum(uint8_t** codeword, int64_t prime, int64_t n, int64_t redundancy, int64_t packet_size);

/// @brief 基于典型值的解码算法
/// @param codeword 码字
/// @param prime 一个码元包含的符号个数
/// @param n 码元个数
/// @param redundancy 冗余度 
/// @param packet_size 符号大小
/// @param erasure_index 失效列下标数组
/// @param len 失效列数量，要求不超过冗余度
void syndrome_based_decode(uint8_t** codeword, int64_t prime, int64_t n, int64_t redundancy, int64_t packet_size, int64_t* erasure_index, int64_t len);

/// @brief 优化的基于典型值的解码算法
/// @param codeword 码字
/// @param prime 一个码元包含的符号个数
/// @param n 码元数量
/// @param redundancy 冗余度
/// @param packet_size 符号大小
/// @param erasure_index 失效列下标数组
/// @param len 失效列数量，要求不超过冗余度
void modified_syndrome_based_decode(uint8_t** codeword, int64_t prime, int64_t n, int64_t redundancy, int64_t packet_size, int64_t* erasure_index, int64_t len);

/// @brief 基于多项式插值的解码算法
/// @param codeword 码字
/// @param prime 一个码元包含的符号个数
/// @param n 码元数量
/// @param redundancy 冗余度
/// @param packet_size 符号大小
/// @param erasure_index 失效列下标数组
/// @param len 失效列数量，要求不超过冗余度
void interpolation_based_decode(uint8_t** codeword, int64_t prime, int64_t n, int64_t redundancy, int64_t packet_size, int64_t* erasure_index, int64_t len);

/// @brief 优化的基于多项式插值的解码算法
/// @param codeword 码字
/// @param prime 一个码元包含的符号个数
/// @param n 码元数量
/// @param redundancy 冗余度
/// @param packet_size 符号大小
/// @param erasure_index 失效列下标数组
/// @param len 失效列数量，要求不超过冗余度
void modified_interpolation_based_decode(uint8_t** codeword, int64_t prime, int64_t n, int64_t redundancy, int64_t packet_size, int64_t* erasure_index, int64_t len);

/// @brief 基于范德蒙德矩阵LU分解的解码算法
/// @param codeword 码字
/// @param prime 一个码元包含的符号个数
/// @param n 码元数量
/// @param redundancy 冗余度
/// @param packet_size 符号大小
/// @param erasure_index 失效列下标数组
/// @param len 失效列数量，要求不超过冗余度
void LU_based_decode(uint8_t** codeword, int64_t prime, int64_t n, int64_t redundancy, int64_t packet_size, int64_t* erasure_index, int64_t len);
#endif