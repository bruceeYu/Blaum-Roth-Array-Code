
#include"arraycode.h"

int64_t mod(int64_t value, int64_t prime) {
    int64_t ratio = value / prime;
    if (value < 0) {
        ratio = (value - prime + 1) / prime;
    }
    return (value - ratio * prime);
}

int64_t max(int64_t a, int64_t b) {
    return (a > b ? a : b);
}

int64_t min(int64_t a, int64_t b) {
    return (a < b ? a : b);
}

int64_t diff(int64_t a, int64_t b) {
    return (max(a, b) - min(a, b));
}

void swap(int64_t* pa, int64_t* pb) {
    int64_t temp = *pa;
    *pa = *pb;
    *pb = temp;
}

uint8_t* locate_packet(const uint8_t* ptr, int64_t index, int64_t prime, int64_t packet_size) {
    index = mod(index, prime);
    return (uint8_t*)ptr + index * packet_size;
}

void extract_commom_factor(int64_t* pairs, const int64_t len, const int64_t prime) {
    int64_t term_exp = 0;
    for (int64_t i = 1; i < len; ++i) {
        term_exp += min(pairs[0], pairs[i]);
        pairs[i] = diff(pairs[0], pairs[i]);
    }
    pairs[0] = term_exp;
}

int64_t integrate(int64_t* pairs, const int64_t len, const int64_t prime) {
    // 以 mod 1+x^7为例，对于1+x 1+x^2 1+x^3 1+x^4=x^4(1+x^3) 1+x^5=x^5(1+x^2) 1+x^6=x^6(1+x) 只需记录1+x 1+x^2 1+x^3 
    const int64_t upper_bound = (prime - 1) / 2;
    int64_t* freq = (int64_t*)calloc(1 + upper_bound, sizeof(int64_t));
    // freq[0]记录单项式x的次数，freq[1]记录1+x的次数，freq[2]记录1+x^2的次数···

    for (int64_t i = 1; i < len; ++i) {
        freq[0] += min(pairs[0], pairs[i]);
        pairs[i] = diff(pairs[0], pairs[i]);

        const int64_t ell = pairs[i];
        if (ell <= upper_bound) {
            freq[ell]++;
        }
        else { // 翻转 1+x^ell = x^ell(1+x^(prime-ell)) mod 1+x^prime
            freq[0] += ell;
            freq[prime - ell]++;
        }
    }

    while (true) {
        for (int64_t ell = 1; ell <= upper_bound; ++ell) {
            if (freq[ell] <= 1) {
                continue;
            }
            else {
                int64_t next = 2 * ell;
                if (next <= upper_bound) {
                    freq[next] += freq[ell] / 2;
                }
                else { // 翻转 1+x^next = x^next(1+x^(prime-next)) mod 1+x^prime
                    freq[0] += next * (freq[ell] / 2);
                    freq[prime - next] += (freq[ell] / 2);
                }
                // 偶数个因式完全合并，奇数个因式合并后剩一项，用位运算实现该逻辑
                freq[ell] &= 1;
            }
        }
        // 检查1+x 1+x^2 ···出现次数是否均不超过1次
        int cnt = 0;
        for (int ell = 1; ell <= upper_bound; ++ell) {
            if (freq[ell] <= 1) {
                ++cnt;
            }
        }
        if (cnt == upper_bound) {
            break;
        }
    }

    int nonzero = 1;
    // x^freq[0] = x^(freq[0] mod prime)
    pairs[0] = mod(freq[0], prime);
    for (int ell = 1; ell <= upper_bound; ++ell) {
        if (freq[ell] > 0) {
            pairs[nonzero++] = ell;
        }
    }
    free(freq);
    freq = NULL;
    return nonzero; // pairs所对应的因式个数 x^pairs[0](1+x^pairs[1])(1+x^pairs[2])...
}

uint8_t** create_codeword(int64_t prime, int64_t n, int64_t redundancy, int64_t packet_size) {

    uint8_t** codeword = (uint8_t**)malloc(n * sizeof(uint8_t*));
    int64_t k = n - redundancy;
    int64_t chunk_size = prime * packet_size;
    for (int64_t j = 0; j < k; ++j) {
        posix_memalign((void**)&codeword[j], 32, chunk_size);
        rand_chunk(codeword[j], prime, packet_size);
    }
    for (int64_t j = k; j < n; ++j) {
        posix_memalign((void**)&codeword[j], 32, chunk_size);
        memset(codeword[j], 0, chunk_size);
    }
    return codeword;
}

uint8_t** destroy_codeword(uint8_t** codeword, int64_t n) {
    for (int64_t j = 0; j < n; ++j) {
        free(codeword[j]);
        codeword[j] = NULL;
    }
    return NULL;
}

uint8_t** create_elements(int64_t amount, int64_t chunk_size) {

    uint8_t** ptr = (uint8_t**)malloc(amount * sizeof(uint8_t*));
    for (int64_t j = 0; j < amount; ++j) {
        posix_memalign((void**)&ptr[j], 32, chunk_size);
        memset((void*)ptr[j], 0, chunk_size);
    }
    return ptr;
}

uint8_t** destroy_elements(uint8_t** mem, int64_t amount) {
    for (int64_t j = 0; j < amount; ++j) {
        free(mem[j]);
        mem[j] = NULL;
    }
    free(mem);
    mem = NULL;
    return NULL;
}

void init_chunk(uint8_t* ptr, int64_t chunk_size) {
    memset(ptr, 0, chunk_size);
}

void rand_chunk(uint8_t* ptr, int64_t prime, int64_t packet_size) {

    for (int64_t i = 0; i < (prime - 1) * packet_size; ++i) {
        *(ptr + i) = rand() % 256;
    }
    memset(ptr + (prime - 1) * packet_size, 0, packet_size);
}

void rectify(uint8_t* ptr, int64_t prime, int64_t packet_size) {
    uint8_t* last = locate_packet(ptr, prime - 1, prime, packet_size);
    if (equal_zero(last, packet_size)) {
        return;
    }

    // symbol_0 -> symbol_{prime - 2}
    for (int64_t i = 0; i < prime - 1; ++i) {

        // symbol_i = symbol_i ^ symbol_{prime-1}
        // 或者 symbol_i ^= symbol_{prime-1}

        uint8_t* arr[3] = { NULL };
        uint8_t* temp = NULL;
        posix_memalign((void**)&temp, 32, packet_size);
        arr[2] = temp;

        arr[0] = locate_packet(ptr, i, prime, packet_size);
        arr[1] = last;
        xor_gen(3, packet_size, (void**)arr);
        memcpy(arr[0], temp, packet_size); // write back
        free(temp);
        temp = NULL;
    }

    // mod Mpx
    memset(last, 0, packet_size);
}

void local_cyclic_shift(uint8_t* ptr, int64_t shift, int64_t prime, int64_t packet_size) {

    shift = mod(shift, prime);
    if (0 == shift) {
        return;
    }

    // 暂存第0个包
    uint8_t* temp = (uint8_t*)malloc(packet_size);
    memcpy(temp, ptr, packet_size);

    int64_t i = 0;
    while (mod(i - shift, prime) != 0) {
        // 第i个包 ⬅ 第prev个包
        int64_t prev = mod(i - shift, prime);
        // overwrite forward!
        memcpy(ptr + i * packet_size, ptr + prev * packet_size, packet_size);
        // iterate backward!
        i = prev;
    }
    memcpy(ptr + i * packet_size, temp, packet_size);
    free(temp);
    temp = NULL;
}

void add_equal(uint8_t* dest, uint8_t* src, int64_t shift, int64_t prime, int64_t packet_size) {
    uint8_t* temp = NULL;
    posix_memalign((void**)&temp, 32, prime * packet_size);
    vector_add(temp, dest, 0, src, shift, prime, packet_size);
    memcpy(dest, temp, prime * packet_size);
    free(temp);
    temp = NULL;
}

void divide_equal(uint8_t* dest, int64_t exp_first, int64_t exp_second, int64_t prime, int64_t packet_size, bool fast) {
    uint8_t* temp = NULL;
    posix_memalign((void**)&temp, 32, prime * packet_size);
    indirect_divide(temp, dest, prime, packet_size, exp_first, exp_second, fast);
    memcpy(dest, temp, prime * packet_size);
    free(temp);
    temp = NULL;
}

void vector_add(uint8_t* dest, const uint8_t* src1, int64_t shift1, const uint8_t* src2, int64_t shift2, int64_t prime, int64_t packet_size) {
    uint8_t* arr[3] = { NULL };
    for (int64_t i = 0; i < prime; ++i) { // ith symbol
        arr[2] = locate_packet(dest, i, prime, packet_size);

        arr[0] = locate_packet(src1, i - shift1, prime, packet_size);
        arr[1] = locate_packet(src2, i - shift2, prime, packet_size);
        xor_gen(3, packet_size, (void**)arr);
    }
}

void indirect_series_mul(uint8_t* gx, const uint8_t* fx, int64_t prime, int64_t packet_size, int64_t* pairs, int64_t len) {

    const uint64_t chunk_size = prime * packet_size;
    memcpy(gx, fx, chunk_size);

    // 顺序迭代 gx = (1 + x^pairs[s])gx = gx + x^pairs[s] gx
    for (int64_t s = 1; s < len; ++s) { // 第s个因式乘法 1+x^pairs[s]
        uint8_t* temp = NULL;
        posix_memalign((void**)&temp, 32, chunk_size);
        for (int64_t i = 0; i < prime; ++i) { // ith symbol
            uint8_t* arr[3] = { NULL };
            arr[2] = locate_packet(temp, i, prime, packet_size);

            arr[0] = locate_packet(gx, i, prime, packet_size);
            arr[1] = locate_packet(gx, i - pairs[s], prime, packet_size);
            xor_gen(3, packet_size, (void**)arr);
        }
        memcpy(gx, temp, chunk_size); // 将临时结果回写到多项式gx中
        free(temp);
        temp = NULL;
    }

    // 边界处理 gx = x^pairs[0] gx
    local_cyclic_shift(gx, pairs[0], prime, packet_size);
}

void rough_divide(uint8_t* gx, const uint8_t* fx, int64_t prime, int64_t packet_size, int64_t exp_first, int64_t exp_second) {

    // gx = fx / (x^exp_first + x^exp_second)
    int64_t d = diff(exp_first, exp_second);
    // 边界处理3个系数
    memset(locate_packet(gx, prime - 1, prime, packet_size), 0, packet_size);
    memcpy(locate_packet(gx, prime - d - 1, prime, packet_size), locate_packet(fx, prime - 1, prime, packet_size), packet_size);
    memcpy(locate_packet(gx, d - 1, prime, packet_size), locate_packet(fx, d - 1, prime, packet_size), packet_size);
    // 其他 prime-3 个系数通过迭代确定
    for (int64_t i = 1; i <= prime - 3; ++i) { // ell_th symbol
        uint8_t* arr[3] = { NULL };
        arr[2] = locate_packet(gx, prime - (i + 1) * d - 1, prime, packet_size);

        arr[0] = locate_packet(gx, prime - i * d - 1, prime, packet_size);
        arr[1] = locate_packet(fx, prime - i * d - 1, prime, packet_size);
        xor_gen(3, packet_size, (void**)arr);
    }
    local_cyclic_shift(gx, -min(exp_first, exp_second), prime, packet_size);
}

void exact_divide(uint8_t* gx, const uint8_t* fx, int64_t prime, int64_t packet_size, int64_t exp_first, int64_t exp_second) {

    int64_t d = diff(exp_first, exp_second);

    // 边界处理1个系数
    uint8_t** arr = (uint8_t**)malloc(((prime - 1) / 2 + 1) * sizeof(uint8_t*));
    arr[(prime - 1) / 2] = gx;

    for (int64_t i = 2; i <= prime - 1; i += 2) {
        arr[i / 2 - 1] = locate_packet(fx, i * d, prime, packet_size);
    }
    xor_gen((prime - 1) / 2 + 1, packet_size, (void**)arr);
    free(arr);
    arr = NULL;

    // 其他prime-1个系数通过迭代确定
    for (int64_t i = 1; i <= prime - 1; ++i) {
        uint8_t* arr[3];
        arr[2] = locate_packet(gx, i * d, prime, packet_size);

        arr[0] = locate_packet(gx, (i - 1) * d, prime, packet_size);
        arr[1] = locate_packet(fx, i * d, prime, packet_size);
        xor_gen(3, packet_size, (void**)arr);
    }

    local_cyclic_shift(gx, -min(exp_first, exp_second), prime, packet_size);
}

void indirect_divide(uint8_t* gx, const uint8_t* fx, int64_t prime, int64_t packet_size, int64_t exp_first, int64_t exp_second, bool fast) {
    if (fast) {
        rough_divide(gx, fx, prime, packet_size, exp_first, exp_second);
    }
    else {
        exact_divide(gx, fx, prime, packet_size, exp_first, exp_second);
    }
}

void indirect_series_divide(uint8_t* fx, int64_t prime, int64_t packet_size, int64_t* pairs, int64_t len) {

    // 边界的处理，除上x^{pairs[0]}相当于乘以x^{-pairs[0]}
    local_cyclic_shift(fx, -pairs[0], prime, packet_size);

    uint8_t* gx = NULL;
    const int64_t chunk_size = prime * packet_size;
    posix_memalign((void**)&gx, 32, chunk_size);
    // 在连续除法时，必须使用保证结果多项式含有偶数个非零项的除法
    for (int64_t i = 1; i < len - 1; ++i) {
        indirect_divide(gx, fx, prime, packet_size, 0, pairs[i], false);
        memcpy(fx, gx, chunk_size); // write back
    }
    // 最后进行的除法可以使用不保证结果多项式含有偶数个非零项的除法，加快速度
    indirect_divide(gx, fx, prime, packet_size, 0, pairs[len - 1], true);
    memcpy(fx, gx, chunk_size); // write back
    free(gx);
    gx = NULL;
}

void direct_divide(uint8_t* gx, const uint8_t* fx, int64_t prime, int64_t packet_size, int64_t exp_first, int64_t exp_second) {
    // set virtual coefficient zero
    memset(locate_packet(gx, prime - 1, prime, packet_size), 0, packet_size);

    uint8_t* arr[4] = { NULL };
    int64_t d = exp_first - exp_second;
    // total prime-1 coefficients
    for (int64_t s = 1; s < prime; ++s) {
        arr[3] = locate_packet(gx, -2 * s * d - 1, prime, packet_size); // dest

        arr[0] = locate_packet(gx, -(2 * s - 2) * d - 1, prime, packet_size); // src1
        arr[1] = locate_packet(fx, -(2 * s - 2) * d + exp_second - 1, prime, packet_size); // src2
        arr[2] = locate_packet(fx, -(2 * s - 1) * d + exp_second - 1, prime, packet_size); // src3
        xor_gen(4, packet_size, (void**)arr);
    }
}

void direct_series_divide(uint8_t* fx, int64_t prime, int64_t packet_size, const int64_t* pairs, int64_t len) {

    const int64_t chunk_size = prime * packet_size;
    uint8_t* gx = NULL;
    posix_memalign((void**)&gx, 32, chunk_size);
    for (int64_t i = 1; i < len; ++i) {
        direct_divide(gx, fx, prime, packet_size, pairs[0], pairs[i]);
        memcpy(fx, gx, chunk_size); // write back
    }
    free(gx);
    gx = NULL;
}

void fill_survival_index(const int64_t* erasure_index, int64_t lambda, int64_t* survival_index, int64_t delta) {

    int64_t* map = (int64_t*)calloc(lambda + delta, sizeof(int64_t));
    for (int64_t i = 0; i < lambda; ++i) {
        int64_t ei = erasure_index[i];
        map[ei] = 1;
    }

    int64_t k = 0;
    int64_t n = lambda + delta;
    for (int64_t j = 0; j < n; ++j) {
        if (!map[j]) {
            survival_index[k++] = j;
        }
    }
    free(map);
    map = NULL;
}

int64_t cal_combination(int64_t n, int64_t m) {
    if (m == 0 || m == n) { // 递归出口 C(n, 0) = C(n, n) = 1
        return 1;
    }
    if (m > n / 2) { // 组合数的互补性质，例如 C(7, 5) = C(7, 2)
        m = n - m; // 减少递归深度 
    }
    // 组合数的定义 C(n, m) = n * C(n - 1, m - 1) / m
    return n * cal_combination(n - 1, m - 1) / m;
}

uint8_t** cal_syndrome(uint8_t** codeword, int64_t prime, int64_t n, int64_t amount, int64_t packet_size, int64_t* erasure_index, int64_t len) {

    const int64_t chunk_size = prime * packet_size;
    uint8_t** syn = (uint8_t**)malloc(amount * sizeof(uint8_t*));
    for (int64_t ell = 0; ell < amount; ++ell) {
        posix_memalign((void**)&syn[ell], 32, chunk_size);
        memset(syn[ell], 0, chunk_size);
    }

    int64_t delta = n - len;
    int64_t* survival_index = (int64_t*)malloc(delta * sizeof(int64_t));
    fill_survival_index(erasure_index, len, survival_index, delta);

    // 计算第ell个典型值，共amount个，要求amount不超过redundancy
    for (int64_t ell = 0; ell < amount; ++ell) {
        // 确定当前典型值的第i个系数，共prime个
        for (int64_t i = 0; i < prime; ++i) {
            uint8_t** arr = (uint8_t**)malloc((delta + 1) * sizeof(uint8_t*));
            arr[delta] = locate_packet(syn[ell], i, prime, packet_size);

            for (int64_t j = 0; j < delta; ++j) {
                int64_t sj = survival_index[j];
                arr[j] = locate_packet(codeword[sj], i - ell * sj, prime, packet_size);
            }
            // 计算当前斜线的校验和得到第ell个典型值的第i个系数
            xor_gen(delta + 1, packet_size, (void**)arr);
            free(arr);
            arr = NULL;
        }
    }
    free(survival_index);
    survival_index = NULL;
    return syn;
}

void simulate_aux(int64_t begin, int64_t end, int64_t amount, int64_t remain, int64_t* output, int64_t** patterns, int64_t* seq) {

    // 当前范围无需再选，将构造的结果output输出到最终结果pattern中
    if (0 == remain) {
        memcpy(patterns[*seq], output, amount * sizeof(int64_t));
        ++(*seq);
        return;
    }

    // 在范围[i,end]中，需取出remain个数，并且将所有组合输出到patterns中
    for (int64_t i = begin; i <= end - remain + 1; ++i) {
        output[amount - remain] = i; // 1.当前位置取出1个
        simulate_aux(i + 1, end, amount, remain - 1, output, patterns, seq); // 2.后续范围取出remain-1个
    }
}

int64_t simulate(int64_t begin, int64_t end, int64_t amount, int64_t** patterns) {

    int64_t* output = (int64_t*)malloc(amount * sizeof(int64_t));
    int64_t seq = 0;
    simulate_aux(begin, end, amount, amount, output, patterns, (int64_t*)&seq);
    free(output);
    output = NULL;
}

bool equal_zero(uint8_t* ptr, int64_t mem_size) {
    for (int64_t order = 0; order < mem_size; ++order) {
        if (ptr[order]) {
            return false;
        }
    }
    return true;
}

bool checksum(uint8_t** codeword, int64_t prime, int64_t n, int64_t redundancy, int64_t packet_size) {

    bool correct = true;

    uint8_t* temp = NULL;
    posix_memalign((void**)&temp, 32, packet_size);

    // 计算斜率为ell的斜线校验和，共prime条
    for (int64_t ell = 0; (ell < redundancy) && correct; ++ell) {

        // 计算斜率为ell的第i条斜线的校验和
        for (int64_t i = 0; i < prime; ++i) {

            uint8_t** arr = (uint8_t**)malloc((n + 1) * sizeof(uint8_t*));
            arr[n] = temp;

            for (int64_t j = 0; j < n; ++j) {
                arr[j] = locate_packet(codeword[j], i - ell * j, prime, packet_size);
            }
            xor_gen(n + 1, packet_size, (void**)arr);
            free(arr);
            arr = NULL;
            // 验证当前斜线的校验和是否为零
            if (!equal_zero(temp, packet_size)) {
                correct = false;
                break;
            }
        }
    }
    free(temp);
    temp = NULL;
    return correct;
}

void syndrome_based_decode(uint8_t** codeword, int64_t prime, int64_t n, int64_t redundancy, int64_t packet_size, int64_t* erasure_index, int64_t len) {

    const int64_t chunk_size = prime * packet_size;

    // 第一部分 计算典型值
    uint8_t** syndrome = cal_syndrome(codeword, prime, n, redundancy, packet_size, erasure_index, len);

    // 第二部分 计算乘积
    uint8_t** Qx = (uint8_t**)malloc((redundancy + len) * sizeof(uint8_t*));
    for (int64_t ell = 0; ell < redundancy; ++ell) {
        Qx[ell] = syndrome[ell];
        syndrome[ell] = NULL;
    }
    free(syndrome);
    syndrome = NULL;
    for (int64_t ell = redundancy; ell < redundancy + len; ++ell) {
        posix_memalign((void**)&Qx[ell], 32, chunk_size);
    }

    for (int64_t s = 0; s < len; ++s) {
        // 特殊处理边界 Qx[redundancy+s] = x^erasure_index[s] Qx[redundancy + s - 1]
        for (int i = 0; i < prime; ++i) {
            uint8_t* dest = locate_packet(Qx[redundancy + s], i, prime, packet_size);
            uint8_t* src = locate_packet(Qx[redundancy + s - 1], i - erasure_index[s], prime, packet_size);
            memcpy(dest, src, packet_size);
        }

        // 逆序迭代 Qx[ell] = Qx[ell] + x^erasure_index[s] Qx[ell-1] 或 Qx[ell] += x^erasure_index[s] Qx[ell-1]
        for (int64_t ell = redundancy - 1 + s; ell >= 1; --ell) {
            add_equal(Qx[ell], Qx[ell - 1], erasure_index[s], prime, packet_size);
        }
    }

    // 第三部分 计算sigma[i]
    uint8_t** sigma = (uint8_t**)malloc(len * sizeof(uint8_t*));
    for (int64_t ell = 0; ell < len; ++ell) {
        posix_memalign((void**)&sigma[ell], 32, chunk_size);
        memcpy(sigma[ell], Qx[0], chunk_size);
    }
    for (int64_t i = 0; i < len; ++i) {
        // sigma[i] = x^erasure_index[i] sigma[i] + Qx[ell]
        for (int64_t ell = 1; ell < len; ++ell) {
            uint8_t* temp = NULL;
            posix_memalign((void**)&temp, 32, chunk_size);
            vector_add(temp, sigma[i], erasure_index[i], Qx[ell], 0, prime, packet_size);
            memcpy(sigma[i], temp, chunk_size);
            free(temp);
            temp = NULL;
        }
        // sigma[i] = sigma[i] mod M_px
        rectify(sigma[i], prime, packet_size);
    }
    Qx = destroy_elements(Qx, redundancy + len);

    // 第四部分 计算除法
    for (int64_t i = 0; i < len; ++i) {
        int64_t ei = erasure_index[i];
        memcpy(codeword[ei], sigma[i], chunk_size);
    }
    sigma = destroy_elements(sigma, len);

    for (int64_t i = 0; i < len; ++i) {
        const int64_t curr_repair = erasure_index[i];
        swap(&erasure_index[0], &erasure_index[i]);
        direct_series_divide(codeword[curr_repair], prime, packet_size, erasure_index, len);
        swap(&erasure_index[0], &erasure_index[i]);
    }
}

void modified_syndrome_based_decode(uint8_t** codeword, int64_t prime, int64_t n, int64_t redundancy, int64_t packet_size, int64_t* erasure_index, int64_t len) {

    const int64_t chunk_size = prime * packet_size;

    // 第一部分 计算 lambda 典型值
    uint8_t** syndrome = cal_syndrome(codeword, prime, n, len, packet_size, erasure_index, len);

    // 第二部分 计算乘积
    uint8_t** Qx = (uint8_t**)malloc(len * sizeof(uint8_t*));
    for (int64_t ell = 0; ell < len; ++ell) {
        Qx[ell] = syndrome[ell];
        syndrome[ell] = NULL;
    }
    free(syndrome);
    syndrome = NULL;

    for (int64_t s = 0; s < len; ++s) {
        // 逆序迭代 Qx[ell] = Qx[ell] + x^erasure_index[s] Qx[ell-1] 或 Qx[ell] += x^erasure_index[s] Qx[ell-1]
        for (int64_t ell = len - 1; ell >= 1; --ell) {
            add_equal(Qx[ell], Qx[ell - 1], erasure_index[s], prime, packet_size);
        }
    }

    // 第三部分 计算sigma
    uint8_t** sigma = create_elements(len, chunk_size);
    for (int64_t i = 0; i < len;++i) {
        memcpy(sigma[i], Qx[0], chunk_size);
    }
    for (int64_t i = 0; i < len; ++i) {
        // sigma[i] = x^erasure_index[i] sigma[i] + Qx[ell]
        for (int64_t ell = 1; ell < len; ++ell) {
            uint8_t* temp = NULL;
            posix_memalign((void**)&temp, 32, chunk_size);
            vector_add(temp, sigma[i], erasure_index[i], Qx[ell], 0, prime, packet_size);
            memcpy(sigma[i], temp, chunk_size);
            free(temp);
            temp = NULL;
        }
    }

    Qx = destroy_elements(Qx, len);

    // 第四部分 计算除法
    for (int64_t i = 0; i < len; ++i) {
        int64_t ei = erasure_index[i];
        memcpy(codeword[ei], sigma[i], chunk_size);
    }
    sigma = destroy_elements(sigma, len);

    for (int64_t i = 0; i < len; ++i) {
        int64_t ei = erasure_index[i];
        int64_t* pairs = (int64_t*)malloc(len * sizeof(int64_t));
        swap(&erasure_index[0], &erasure_index[i]);
        memcpy(pairs, erasure_index, len * sizeof(int64_t));
        int64_t nonzero = integrate(pairs, len, prime);
        indirect_series_divide(codeword[ei], prime, packet_size, pairs, nonzero);
        swap(&erasure_index[0], &erasure_index[i]);
        free(pairs);
        pairs = NULL;
    }
}

void interpolation_based_decode(uint8_t** codeword, int64_t prime, int64_t n, int64_t redundancy, int64_t packet_size, int64_t* erasure_index, int64_t len) {

    const int64_t chunk_size = prime * packet_size;
    int64_t delta = n - len;
    int64_t* survival_index = (int64_t*)malloc(delta * sizeof(int64_t));
    fill_survival_index(erasure_index, len, survival_index, delta);

    // 第一部分，计算delta个a_ploys[j]
    uint8_t** a_ploys = create_elements(delta, chunk_size);

    for (int64_t j = 0; j < delta; ++j) {
        // a_ploys[j] = (x^hj + x^erasure_index[0])(x^hj + erasure_index[1])...(x^hj + x^erasure_index[len-1]) codeword[hj]
        const int64_t hj = survival_index[j];
        // 因式连乘相当于 (x^pairs[0] + x^pairs[1])(x^pairs[0] + x^pairs[1])...(x^pairs[0] + x^pairs[len])
        int64_t* pairs = (int64_t*)malloc((len + 1) * sizeof(int64_t));
        pairs[0] = hj;
        memcpy(pairs + 1, erasure_index, len * sizeof(int64_t)); // factor (x^hj + x^ei)

        // 对于每个因式提取公共项，相当于x^pairs[0](1 + x^pairs[1])(1 + x^pairs[2])...(1 + x^pairs[len])
        extract_commom_factor(pairs, len + 1, prime);
        indirect_series_mul(a_ploys[j], codeword[hj], prime, packet_size, pairs, len + 1);
        rectify(a_ploys[j], prime, packet_size); // mod Mpx！！！
        free(pairs);
        pairs = NULL;
    }

    // 第2部分，计算len个b_ploys[i]
    uint8_t** b_ploys = create_elements(len, chunk_size);

    for (int i = 0; i < len; ++i) {
        // b_ploys[i] = a_ploys[0] / (x^survival_index[0] + x^erasure_index[i]) + ... + a_ploys[delta-1] / (x^survival_index[delta-1] + x^erasure_index[i])
        const int64_t ei = erasure_index[i];

        uint8_t** arr = (uint8_t**)malloc((delta + 1) * sizeof(uint8_t*));
        arr[delta] = b_ploys[i];

        uint8_t** temp_ploys = create_elements(delta, chunk_size);
        for (int j = 0; j < delta; ++j) {
            const int64_t hj = survival_index[j];
            direct_divide(temp_ploys[j], a_ploys[j], prime, packet_size, hj, ei);
            arr[j] = temp_ploys[j];
        }
        xor_gen(delta + 1, chunk_size, (void**)arr);

        free(arr);
        arr = NULL;

        temp_ploys = destroy_elements(temp_ploys, delta);
    }

    free(survival_index);
    survival_index = NULL;

    a_ploys = destroy_elements(a_ploys, delta);

    // 第3部分 计算连续除法
    for (int64_t i = 0; i < len; ++i) {
        int64_t ei = erasure_index[i];
        memcpy(codeword[ei], b_ploys[i], chunk_size);
    }
    b_ploys = destroy_elements(b_ploys, len);

    for (int64_t i = 0; i < len; ++i) {
        const int64_t curr_repair = erasure_index[i];
        swap(&erasure_index[0], &erasure_index[i]);
        direct_series_divide(codeword[curr_repair], prime, packet_size, erasure_index, len);
        swap(&erasure_index[0], &erasure_index[i]);
    }
}

void modified_interpolation_based_decode(uint8_t** codeword, int64_t prime, int64_t n, int64_t redundancy, int64_t packet_size, int64_t* erasure_index, int64_t len) {

    const int64_t chunk_size = prime * packet_size;
    int64_t delta = n - len;
    int64_t* survival_index = (int64_t*)malloc(delta * sizeof(int64_t));
    fill_survival_index(erasure_index, len, survival_index, delta);

    // 第一部分，计算delta个a_ploys[j] 
    uint8_t** a_ploys = create_elements(delta, chunk_size);

    for (int64_t j = 0; j < delta; ++j) {
        // a_ploys[j] = (x^hj + x^erasure_index[0])(x^hj + erasure_index[1])...(x^hj + x^erasure_index[len-1]) codeword[hj]
        const int64_t hj = survival_index[j];
        // 因式连乘相当于 (x^pairs[0] + x^pairs[1])(x^pairs[0] + x^pairs[1])...(x^pairs[0] + x^pairs[len])
        int64_t* pairs = (int64_t*)malloc((len + 1) * sizeof(int64_t));
        pairs[0] = hj;
        memcpy(pairs + 1, erasure_index, len * sizeof(int64_t)); // factor (x^hj + x^ei)

        // 整理合并因式连乘，相当于x^pairs[0](1 + x^pairs[1])(1 + x^pairs[2])...(1 + x^pairs[nonzero-1])
        int64_t nonzero = integrate(pairs, len + 1, prime);

        indirect_series_mul(a_ploys[j], codeword[hj], prime, packet_size, pairs, nonzero);
        free(pairs);
        pairs = NULL;
    }

    // 第2部分，计算len个b_ploys[i]
    uint8_t** b_ploys = create_elements(len, chunk_size);

    for (int i = 0; i < len; ++i) { // ith b_i(x)
        // b_ploys[i] = a_ploys[0] / (x^survival_index[0] + x^erasure_index[i]) + ... + a_ploys[delta-1] / (x^survival_index[delta-1] + x^erasure_index[i])
        uint8_t** arr = (uint8_t**)malloc((delta + 1) * sizeof(uint8_t*));
        arr[delta] = b_ploys[i];

        uint8_t** temp_ploys = create_elements(delta, chunk_size);
        for (int j = 0; j < delta; ++j) { // a_ploys[j] / (x^survival_index[j] + x^erasure_index[i])
            exact_divide(temp_ploys[j], a_ploys[j], prime, packet_size, survival_index[j], erasure_index[i]);
            arr[j] = temp_ploys[j];
        }
        xor_gen(delta + 1, chunk_size, (void**)arr);

        free(arr);
        arr = NULL;
        temp_ploys = destroy_elements(temp_ploys, delta);
    }

    free(survival_index);
    survival_index = NULL;

    a_ploys = destroy_elements(a_ploys, delta);

    // 第3部分 计算连续除法
    for (int64_t i = 0; i < len; ++i) {
        memcpy(codeword[erasure_index[i]], b_ploys[i], chunk_size);
    }
    b_ploys = destroy_elements(b_ploys, len);

    for (int64_t i = 0; i < len; ++i) {
        int64_t ei = erasure_index[i];
        int64_t* pairs = (int64_t*)malloc(len * sizeof(int64_t));
        swap(&erasure_index[0], &erasure_index[i]);
        memcpy(pairs, erasure_index, len * sizeof(int64_t));
        int64_t nonzero = integrate(pairs, len, prime);
        indirect_series_divide(codeword[ei], prime, packet_size, pairs, nonzero);
        swap(&erasure_index[0], &erasure_index[i]);
        free(pairs);
        pairs = NULL;
    }
}

void LU_based_decode(uint8_t** codeword, int64_t prime, int64_t n, int64_t redundancy, int64_t packet_size, int64_t* erasure_index, int64_t len) {
    const int64_t chunk_size = prime * packet_size;
    // 计算典型值 
    uint8_t** syndrome = cal_syndrome(codeword, prime, n, len, packet_size, erasure_index, len);
    for (int64_t i = 0; i < len; ++i) {
        memcpy(codeword[erasure_index[i]], syndrome[i], chunk_size);
    }
    syndrome = destroy_elements(syndrome, len);


    // 消去上三角矩阵
    for (int64_t ell = 1; ell < len; ++ell) {
        for (int64_t i = len - ell; i < len;++i) {
            // erasure_index[i] 简写为e_i
            // codeword[e_i] = codeword[e_i] + x^e_{ell+i-len} codeword[e_{i-1}] 或 codeword[e_i] += x^e_{ell+i-len} codeword[e_{i-1}]
            add_equal(codeword[erasure_index[i]], codeword[erasure_index[i - 1]], erasure_index[ell + i - len], prime, packet_size);
        }
    }

    // 消去下三角矩阵
    for (int64_t ell = len - 1; ell >= 1; --ell) {
        // 一、处理右边界

        // erasure_index[i] 简写为 e_i
        // codeword[e_{len-1}] = codeword[e_{len-1}] / (x^e_{len - 1} + x^e_{len - ell - 1}) 
        // 或者 codeword[e_{len - 1}] /= (x^e_{len - 1} + x^e_{len - ell - 1})
        divide_equal(codeword[erasure_index[len - 1]], erasure_index[len - 1], erasure_index[len - ell - 1], prime, packet_size, (1 == ell));

        // 二、处理区间
        for (int64_t i = len - 2; i >= len - ell; --i) {
            // erasure_index[i] 简写为 e_i
            // codeword[e_i] = (codeword[e_i] - codeword[e_{i + 1}]) / (x^e_i + x^e_{len - ell - 1})
            uint8_t* temp_sub;
            posix_memalign((void**)&temp_sub, 32, chunk_size);
            vector_add(temp_sub, codeword[erasure_index[i]], 0, codeword[erasure_index[i + 1]], 0, prime, packet_size);

            indirect_divide(codeword[erasure_index[i]], temp_sub, prime, packet_size, erasure_index[i], erasure_index[len - ell - 1], (i == len - ell));
            free(temp_sub);
            temp_sub = NULL;
        }

        // 三、处理左边界
        // erasure_index[i] 简写为 e_i
        // codeword[e_{len - ell - 1}] = codeword[e_{len - ell - 1}] - codeword[e_{len - ell}]
        // 或者 codeword[e_{len - ell - 1}] -= x^0 codeword[e_{len - ell}]
        add_equal(codeword[erasure_index[len - ell - 1]], codeword[erasure_index[len - ell]], 0, prime, packet_size);

        // 此时完成消去1个下三角矩阵
    }
}
