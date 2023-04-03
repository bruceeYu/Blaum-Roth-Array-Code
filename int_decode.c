
#include"arraycode.h"
#include <time.h>


void encode(uint8_t** codeword, int64_t prime, int64_t n, int64_t redundancy, int64_t packet_size) {

    int64_t* erasure_index = (int64_t*)malloc(redundancy * sizeof(int64_t));
    for (int64_t s = 0; s < redundancy; ++s) {
        erasure_index[s] = n - redundancy + s;
    }
    interpolation_based_decode(codeword, prime, n, redundancy, packet_size, erasure_index, redundancy);
    free(erasure_index);
    erasure_index = NULL;
}

int main() {
    srand((unsigned)time(NULL));
    int64_t choice = 1;
    while (choice) {
        printf("stripes, prime, n, redundancy, erasures, packet_size(MB)\n");
        int64_t stripes;
        int64_t p, n, r;
        int64_t packet_size; // 符号（包）大小，以字节为单位
        int64_t lambda;
        scanf("%ld %ld %ld %ld %ld %ld", &stripes, &p, &n, &r, &lambda, &packet_size);

        packet_size *= (1024 * 1024); // Byte
        int64_t k = n - r;
        int64_t chunk_size = p * packet_size;

        uint8_t** codeword = create_codeword(p, n, r, packet_size);

        clock_t start_encode = clock();
        for (int64_t i = 0; i < stripes; ++i) {

            encode(codeword, p, n, r, packet_size);
            if (checksum(codeword, p, n, r, packet_size)) {
                printf("encode success!\n");
            }
            else {
                printf("encode fail!\n");
            }
        }
        clock_t finish_encode = clock();

        int64_t encode_data_volume = stripes * r * chunk_size / 1024 / 1024;
        double duration_encode = (double)(finish_encode - start_encode) / CLOCKS_PER_SEC;
        double speed = encode_data_volume / duration_encode;
        printf("encode data volumn(MB): %ld\n", encode_data_volume);
        printf("encode time(s): %lf\n", duration_encode);
        printf("encode speed(MB/s): %lf\n\n", speed);


        // decode
        int64_t total = cal_combination(n, lambda);
        int64_t** patterns = (int64_t**)malloc(total * sizeof(int64_t*));
        for (int64_t i = 0; i < total; ++i) {
            patterns[i] = (int64_t*)calloc(lambda, sizeof(int64_t));
        }

        simulate(0, n - 1, lambda, patterns);

        clock_t start_decode = clock();
        for (int64_t s = 0; s < stripes; ++s) {
            for (int64_t i = 0; i < total; ++i) {
                interpolation_based_decode(codeword, p, n, r, packet_size, patterns[i], lambda);
                if (checksum(codeword, p, n, r, packet_size)) {
                    printf("decode success!\n");
                }
                else {
                    printf("decode fail!\n");
                }
            }
        }
        clock_t finish_decode = clock();

        for (int64_t j = 0; j < total; ++j) {
            free(patterns[j]);
            patterns[j] = NULL;
        }
        free(patterns);
        patterns = NULL;

        int64_t decode_data_volume = stripes * total * lambda * chunk_size / 1024 / 1024;
        double duration_decode = (double)(finish_decode - start_decode) / CLOCKS_PER_SEC;
        double decode_speed = decode_data_volume / duration_decode;
        printf("decode data volumn(MB): %ld\n", decode_data_volume);
        printf("decode time(s): %lf\n", duration_decode);
        printf("decode speed(MB/s): %lf\n\n", decode_speed);

        codeword = destroy_codeword(codeword, n);

        printf("End(0) or continue(1):\n");
        scanf("%ld", &choice);
    }
    return 0;
}