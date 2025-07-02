/*
 This class use:
https://github.com/bshoshany/thread-pool
https://github.com/ReneNyffenegger/cpp-base64
 */
#ifndef GOLDWASSER_MICALI_MULTITHREAD_H
#define GOLDWASSER_MICALI_MULTITHREAD_H
#include <map>
#include <string>
#include <vector>


class Goldwasser_Micali_Multithread {
    public:
    static bool is_prime(unsigned long long &number);
    static void is_prime_for_thread(unsigned long long &number, unsigned long long &prime);
    static unsigned long long generate_positive_number_of_length_gt_0(unsigned int length, bool srand_bool=false);
    static unsigned long long generate_prime_of_length(unsigned int length, bool srand_bool=false);
    static void generate_prime_of_length(unsigned int length,unsigned long long &prime_of_length, bool srand_bool=false, int seed = NULL);
    static unsigned long long generate_prime_of_length_multithread(unsigned int length, bool srand_bool=false);
    static int jacobi_symbol(unsigned long long k, unsigned long long n);
    static int legendre_symbol(unsigned long long a, unsigned long long p);
    static void legendre_symbol(short int &result, unsigned long long a, unsigned long long p);
    static std::map<std::string, unsigned long long> GM_generate_key(int length_of_p, int length_of_q, bool srand_bool=true);
    static std::map<std::string, unsigned long long> GM_generate_key_mutlithread(int length_of_p, int length_of_q, bool srand_bool=true);
    static std::vector<unsigned __int8> string_to_bits(std::string &text);
    static std::string bits_to_string(std::vector<unsigned __int8> &bits);
    static unsigned long long encrypt_bit(unsigned __int8 m_i, unsigned long long x, unsigned long long N);
    static void encrypt_bit(unsigned __int8 m_i, unsigned long long x, unsigned long long N, std::vector<unsigned long long> &ciphertext, unsigned long long i);
    static std::vector<unsigned long long> GM_encode(std::string &plaintext, unsigned long long &x, unsigned long long &N, bool srand_bool=true);
    static std::vector<unsigned long long> GM_encode_multithread(std::string &plaintext, unsigned long long &x, unsigned long long &N, bool srand_bool=true, unsigned short int numberOfThreads = 0);
    static std::string GM_decode(std::vector<unsigned long long> ciphertext, unsigned long long &p, unsigned long long &q);
    static std::string GM_decode_mutlithread(std::vector<unsigned long long> ciphertext, unsigned long long &p, unsigned long long &q, unsigned short int numberOfThreads = 0);
    static void decode_value(std::vector<unsigned __int8> &plaintext_bits, unsigned long long ciphervalue, unsigned long long i, unsigned long long p, unsigned long long q);
};



#endif //GOLDWASSER_MICALI_MULTITHREAD_H
