#include "Goldwasser_Micali_Multithread.h"

#include <assert.h>
#include <bitset>
#include <cmath>
#include <ctime>
#include <numeric>

#include "base64.h"
#include "BS_thread_pool.hpp"

//https://stackoverflow.com/a/14418599
bool Goldwasser_Micali_Multithread::is_prime(unsigned long long &number) {
    if(number < 2) return false;
    if(number == 2) return true;
    if(number % 2 == 0) return false;
    for(unsigned long long i=3; (i*i)<=number; i+=2){
        if(number % i == 0 ) return false;
    }
    return true;
}

void Goldwasser_Micali_Multithread::
is_prime_for_thread(unsigned long long &number, unsigned long long &prime) {
    if(number < 2) return;
    if(number == 2) {
        prime = 2;
        return;
    }
    if(number % 2 == 0) return;
    for(unsigned long long i=3; (i*i)<=number; i+=2){
        if(number % i == 0 || prime != NULL) return;
    }
    prime = number;
}

unsigned long long Goldwasser_Micali_Multithread::generate_positive_number_of_length_gt_0(unsigned int length, bool srand_bool) {
    assert(length > 0 && "Minimum length = 1");
    if (srand_bool) {srand(time({}));}
    unsigned long long random_number;
    random_number = rand() % 10;
    while (random_number == 0) {
        random_number = rand() % 10;
    }
    for (int i = 1; i < length; i++) {
        random_number *= 10;
        random_number += rand() % 10;
    }
    return random_number;
}

unsigned long long Goldwasser_Micali_Multithread::generate_prime_of_length(unsigned int length, bool srand_bool) {
    if (srand_bool) srand(time({}));
    unsigned long long potential_prime;
    potential_prime = generate_positive_number_of_length_gt_0(length);
    if (length == 1) {
        while (is_prime(potential_prime) == false) {
            potential_prime = generate_positive_number_of_length_gt_0(length);
        }
        return potential_prime;
    }
    unsigned long long n;
    unsigned long long n_limit;
    while (true) {
        if (is_prime(potential_prime))
            return potential_prime;
        n = (potential_prime)/2;
        n_limit = pow(10, (length-1));
        if (2*n-1<n_limit) {
            n_limit = 2*n+1;
            if (is_prime(n_limit)) return n_limit;
            n++;
        }
        n_limit = pow(10,length);
        for (n ; 2*n-1<n_limit; n++) {
            potential_prime = 2*n - 1;
            if (is_prime(potential_prime)) return potential_prime;
            potential_prime = 2*n + 1;
            if (potential_prime < n_limit && is_prime(potential_prime)) return potential_prime;
        }
        potential_prime = generate_positive_number_of_length_gt_0(length);
    }
}

void Goldwasser_Micali_Multithread::generate_prime_of_length(unsigned int length,
    unsigned long long &prime_of_length, bool srand_bool, int seed) {
    if (srand_bool) {
        if (seed) srand(seed);
        else srand(time({}));
    }
    prime_of_length = generate_positive_number_of_length_gt_0(length);
    if (length == 1) {
        while (is_prime(prime_of_length) == false) {
            prime_of_length = generate_positive_number_of_length_gt_0(length);
        }
        return;
    }
    unsigned long long n;
    unsigned long long n_limit;
    while (true) {
        if (is_prime(prime_of_length))
            return;
        n = (prime_of_length)/2;
        n_limit = pow(10,(length-1));
        if (2*n-1<n_limit) {
            n_limit = 2*n+1;
            if (is_prime(n_limit)) return;
            n++;
        }
        n_limit = pow(10,length);
        for (n ; 2*n-1<n_limit; n++) {
            prime_of_length = 2*n - 1;
            if (is_prime(prime_of_length)) return;
            prime_of_length = 2*n + 1;
            if (prime_of_length < n_limit && is_prime(prime_of_length)) return;
        }
        prime_of_length = generate_positive_number_of_length_gt_0(length);
    }
}

unsigned long long Goldwasser_Micali_Multithread::generate_prime_of_length_multithread(unsigned int length, bool srand_bool) {
    if (srand_bool) srand(time({}));
    unsigned long long potential_prime;
    potential_prime = generate_positive_number_of_length_gt_0(length);
    if (length == 1) {
        while (is_prime(potential_prime) == false) {
            potential_prime = generate_positive_number_of_length_gt_0(length);
        }
        return potential_prime;
    }
    BS::thread_pool pool((int)(std::thread::hardware_concurrency()*0.8));
    unsigned long long n;
    unsigned long long n_limit;
    unsigned long long prime = NULL;
    // std::future<void> my_future;

    while (true) {
        // my_future = pool.submit_task([&potential_prime, &prime]{is_prime_for_thread(potential_prime, prime);});
        pool.detach_task([&potential_prime, &prime]{is_prime_for_thread(potential_prime, prime);});
        n = (potential_prime)/2;
        n_limit = pow(10,(length-1));
        if (2*n-1<n_limit) {
            n_limit = 2*n+1;
            // my_future = pool.submit_task([&potential_prime, &prime]{is_prime_for_thread(potential_prime, prime);});
            pool.detach_task([&potential_prime, &prime]{is_prime_for_thread(potential_prime, prime);});
            n++;
        }
        n_limit = pow(10,length);
        for (n ; 2*n-1<n_limit && prime != NULL; n++) {
            potential_prime = 2*n - 1;
            // my_future = pool.submit_task([&potential_prime, &prime]{is_prime_for_thread(potential_prime, prime);});
            pool.detach_task([&potential_prime, &prime]{is_prime_for_thread(potential_prime, prime);});
            potential_prime = 2*n + 1;
            if (potential_prime < n_limit) {
                // my_future = pool.submit_task([&potential_prime, &prime]{is_prime_for_thread(potential_prime, prime);});
                pool.detach_task([&potential_prime, &prime]{is_prime_for_thread(potential_prime, prime);});
            }
        }
        while (prime == NULL && pool.get_tasks_total() > 0){}
        if (prime != NULL) {
            pool.purge();
            pool.wait();
            return prime;
        }
        pool.purge();
        pool.wait();
        potential_prime = generate_positive_number_of_length_gt_0(length);
    }
}

//https://rosettacode.org/wiki/Jacobi_symbol#C++
int Goldwasser_Micali_Multithread::jacobi_symbol(unsigned long long k, unsigned long long n) {
    assert(n > 0 && n % 2 == 1 && "n must be odd");
    k %= n;
    int t = 1;
    while (k != 0) {
        while (k % 2 == 0) {
            k /= 2;
            int r = n % 8;
            if (r == 3 || r == 5)
                t = -t;
        }
        std::swap(k, n);
        if (k % 4 == 3 && n % 4 == 3)
            t = -t;
        k %= n;
    }
    return n == 1 ? t : 0;
}

int Goldwasser_Micali_Multithread::legendre_symbol(unsigned long long a, unsigned long long p) {
    assert(is_prime(p));
    assert(p>2);
    if (a == 0) return 0;
    if (a%p == 0) return 0;
    if (a==1) return 1;
    return jacobi_symbol(a,p);
}

void Goldwasser_Micali_Multithread::legendre_symbol(short int &result, unsigned long long a, unsigned long long p) {
    assert(is_prime(p));
    assert(p>2);
    if (a == 0) {
        result = 0;
        return;
    }
    if (a%p == 0) {
        result = 0;
        return;
    }
    if (a==1) {
        result = 1;
        return;
    }
    result = jacobi_symbol(a,p);
}

std::map<std::string, unsigned long long> Goldwasser_Micali_Multithread::GM_generate_key(int length_of_p, int length_of_q, bool srand_bool) {
    if (srand_bool) {
        srand(time({}));
    }
    std::map<std::string, unsigned long long> key;
    key["p"] = generate_prime_of_length(length_of_p);
    key["q"] = generate_prime_of_length(length_of_q);
    while (key["p"]==2 || key["q"]==2 || key["p"]==key["q"]) {
        key["p"] = generate_prime_of_length(length_of_p);
        key["q"] = generate_prime_of_length(length_of_q);
    }
    key["N"] = key["p"]*key["q"];
    //checking overflow
    assert(key["p"]<key["N"] && key["q"]<key["N"] && key["N"]%key["p"]==0 && key["N"]%key["q"]==0);
    unsigned long long x = rand()%key["N"];
    while (x <= 2) {
        x = rand()%key["N"];
    }
    while (legendre_symbol(x,key["p"]) != -1 || legendre_symbol(x,key["q"]) != -1) {
        x = rand()%key["N"];
        while (x <= 2) {
            x = rand()%key["N"];
        }
    }
    assert(jacobi_symbol(x,key["N"]) == 1);
    key["x"] = x;
    return key;
}

std::map<std::string, unsigned long long> Goldwasser_Micali_Multithread::GM_generate_key_mutlithread(int length_of_p,
    int length_of_q, bool srand_bool) {
    if (srand_bool) {
        srand(time({}));
    }
    std::map<std::string, unsigned long long> key;
    unsigned long long * p_a = &key["p"];
    unsigned long long * q_a = &key["q"];
    int seed = rand();
    BS::thread_pool pool((int)(std::thread::hardware_concurrency()*0.8));
    pool.detach_task(
        [length_of_p, p_a, seed]
        {return generate_prime_of_length(length_of_p, *p_a, true, seed);}
    );
    seed = rand();
    pool.detach_task(
            [length_of_q, q_a, seed]
            {return generate_prime_of_length(length_of_q, *q_a, true, seed);}
        );
    pool.wait();
    key["N"] = key["p"]*key["q"];
    assert(key["p"]<key["N"] && key["q"]<key["N"] && key["N"]%key["p"]==0 && key["N"]%key["q"]==0);
    unsigned long long x = rand()%key["N"];
    while (x <= 2) {
        x = rand()%key["N"];
    }
    short int leg1,leg2;
    unsigned long long p = key["p"], q = key["q"];
    while (true) {
        pool.detach_task([&leg1,x, p]
            {return legendre_symbol(leg1,x,p);}
            );
        pool.detach_task([&leg2,x, q]
            {return legendre_symbol(leg2,x,q);}
            );
        pool.wait();
        if (legendre_symbol(x,key["p"]) == -1 && legendre_symbol(x,key["q"]) == -1)
            break;
        x = rand()%key["N"];
        while (x <= 2) {
            x = rand()%key["N"];
        }
    }
    assert(jacobi_symbol(x,key["N"]) == 1);
    key["x"] = x;
    return key;
}

std::vector<unsigned __int8> Goldwasser_Micali_Multithread::string_to_bits(std::string &text) {
    std::string base64_text = base64_encode(text);
    std::vector<unsigned __int8> text_bits;
    std::bitset<8> bits;
    for (int i = 0; i < base64_text.length(); i++) {
        bits = base64_text[i];
        for (__int8 j = size(bits)-1; j >= 0; j--) {
            text_bits.push_back(bits[j]);
        }
    }
    return text_bits;
}

std::string Goldwasser_Micali_Multithread::bits_to_string(std::vector<unsigned char> &bits) {
    std::string base64_text;
    std::bitset<8> tempbitset;
    std::string tempstring;
    for (unsigned long long i = 0; i < bits.size(); i+=8) {
        for (unsigned long long j = i+7; j > i; j--) {
            if (bits[j]) {
                tempbitset.set(7-(j-i));
            }
        }
        base64_text += tempbitset.to_ulong();
        tempbitset.reset();
    }
    return base64_decode(base64_text);
}

unsigned long long Goldwasser_Micali_Multithread::encrypt_bit(unsigned __int8 m_i, unsigned long long x, unsigned long long N) {
    unsigned long long y_i = rand()%N;
    while (std::gcd(y_i,N) != 1) {
        y_i = rand()%N;
    }
    unsigned long long result = (unsigned long long) (pow(y_i,2)*(pow(x,m_i)))%N;
    assert(result<=N);
    assert(pow(y_i,2)>=y_i);
    return result;
}

void Goldwasser_Micali_Multithread::encrypt_bit(unsigned char m_i, unsigned long long x, unsigned long long N,
    std::vector<unsigned long long> &ciphertext, unsigned long long i) {
    unsigned long long y_i = rand()%N;
    while (std::gcd(y_i,N) != 1) {
        y_i = rand()%N;
    }
    unsigned long long result = (unsigned long long) (pow(y_i,2)*(pow(x,m_i)))%N;
    assert(result<=N);
    assert(pow(y_i,2)>=y_i);
    ciphertext[i] = result;
}

std::vector<unsigned long long> Goldwasser_Micali_Multithread::
GM_encode(std::string &plaintext, unsigned long long &x, unsigned long long &N, bool srand_bool) {
    if (srand_bool) {
        srand(time({}));
    }
    std::vector<unsigned __int8> plaintext_bits = Goldwasser_Micali_Multithread::string_to_bits(plaintext);
    std::vector<unsigned long long> ciphertext;
    for (unsigned long long i = 0; i < plaintext_bits.size(); i++) {
        ciphertext.push_back(encrypt_bit(plaintext_bits[i],x,N));
    }
    return ciphertext;
}

std::vector<unsigned long long> Goldwasser_Micali_Multithread::GM_encode_multithread(std::string &plaintext,
    unsigned long long &x, unsigned long long &N, bool srand_bool, unsigned short int numberOfThreads) {
    if (srand_bool) {
        srand(time({}));
    }
    BS::thread_pool pool((int)(std::thread::hardware_concurrency()*0.8));
    if (numberOfThreads != 0) {
        pool.reset(numberOfThreads);
    }
    std::vector<unsigned __int8> plaintext_bits = Goldwasser_Micali_Multithread::string_to_bits(plaintext);
    std::vector<unsigned long long> ciphertext(plaintext_bits.size(),0);
    // BS::thread_pool pool((int)(std::thread::hardware_concurrency()*0.8));
    // for (unsigned long long i = 0; i < plaintext_bits.size(); i++) {
    //     pool.submit_task([plaintext_bits, i, x, N, &ciphertext]
    //         {encrypt_bit(plaintext_bits[i],x,N, ciphertext, i);}
    //         );
    //     //encrypt_bit(plaintext_bits[i],x,N, ciphertext, i);
    // }
    pool.submit_loop(0, plaintext_bits.size(), [plaintext_bits, x, N, &ciphertext](const std::size_t i)
            {encrypt_bit(plaintext_bits[i],x,N, ciphertext, i);}
            );
    pool.wait();
    return ciphertext;
}

std::string Goldwasser_Micali_Multithread::GM_decode(std::vector<unsigned long long> ciphertext, unsigned long long &p,
                                                     unsigned long long &q) {
    std::vector<unsigned __int8> plaintext_bits;
    for (unsigned long long i = 0; i < ciphertext.size(); i++) {
        if (legendre_symbol(ciphertext[i],p) == 1 && legendre_symbol(ciphertext[i],q) == 1) plaintext_bits.push_back(0);
        else plaintext_bits.push_back(1);
    }
    return bits_to_string(plaintext_bits);
}

std::string Goldwasser_Micali_Multithread::GM_decode_mutlithread(std::vector<unsigned long long> ciphertext,
    unsigned long long &p, unsigned long long &q, unsigned short int numberOfThreads) {
    std::vector<unsigned __int8> plaintext_bits(ciphertext.size());
    BS::thread_pool pool((int)(std::thread::hardware_concurrency()*0.8));
    if (numberOfThreads != 0) {
        pool.reset(numberOfThreads);
    }
    pool.submit_loop(0, ciphertext.size(),[&plaintext_bits, ciphertext, p, q] (const std::size_t i)
    {decode_value(plaintext_bits,ciphertext[i], i, p,q);}
    );
    pool.wait();
    return bits_to_string(plaintext_bits);
}

void Goldwasser_Micali_Multithread::decode_value(std::vector<unsigned __int8> &plaintext_bits,
    unsigned long long ciphervalue, unsigned long long i, unsigned long long p, unsigned long long q) {
    if (legendre_symbol(ciphervalue,p) == 1 && legendre_symbol(ciphervalue,q) == 1) {
        plaintext_bits[i] = 0;
    } else plaintext_bits[i] = 1;
}
