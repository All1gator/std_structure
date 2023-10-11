#include <iostream>
#include <cstring>
#include <vector>
#include <string>
#include <cassert>


class BigInteger{
private:
    long long size = 0;
    std::vector<long long> numbers;
    bool sign = 0;
    long long mod = 1e9;

    bool comp_shift(const BigInteger &second, long long shift) {
        for (long long i = second.size + shift; i < size; i++) {
            if (numbers[i] > 0){
                return 1;
            }
        }
        for (long long i = second.size + shift - 1; i >= shift; i--) {
            if (size <= i) {
                return 0;
            }
            if (numbers[i] < second.numbers[i - shift]) {
                return 0;
            }
            if (numbers[i] > second.numbers[i - shift]) {
                return 1;
            }
        }
        return 1;
    }

    void delete_zero() {
        while (size > 1 && numbers[size - 1] == 0) {
            --size;
            numbers.pop_back();
        }
    }

public:
    BigInteger() = default;
    BigInteger(long long number): size(1) {
        sign = ((number < 0) ? 1 : 0);
        if (number < 0) {
            number = -number;
        }
        numbers.push_back(number);
        long long next_digit = numbers[0] / mod;
        numbers[0] %= mod;
        if (next_digit > 0){
            size++;
            numbers.push_back(next_digit);
        }
    }

    void swap(BigInteger& a) {
        std::swap(size, a.size);
        std::swap(sign, a.sign);
        std::swap(numbers, a.numbers);
    }



    BigInteger &operator=(const BigInteger& second) {
        if (&second == this){
            return *this;
        }
        BigInteger copy(second);
        swap(copy);
        return *this;

    }
    bool operator==(const BigInteger& second) const{
        if (size != second.size || (sign != second.sign)){
            return 0;
        }
        for (long long i = 0; i < size; i++){
            if (numbers[i] != second.numbers[i]){
                return 0;
            }
        }
        return 1;
    }

    bool operator!=(const BigInteger& second) const{
        return (!(*this == second));
    }

    bool operator<(const BigInteger& second) const{
        if (sign && second.sign == 0) {
            return 1;
        }
        if (sign == 0 && second.sign) {
            return 0;
        }
        if ((size < second.size && (sign == 0)) || (sign == 1 && size > second.size)) {
            return 1;
        }
        if ((size > second.size && (sign == 0)) || (sign == 1 && size < second.size)) {
            return 0;
        }
        for (long long i = size - 1; i >= 0; i--) {
            if (numbers[i] < second.numbers[i]) {
                return sign ? 0 : 1;
            }
            if (numbers[i] > second.numbers[i]) {
                return sign ? 1 : 0;
            }
        }
        return 0;


    }
    bool operator<=(const BigInteger& second) const {
        return ((*this < second) || (*this == second));

    }

    bool operator>(const BigInteger& second) const {
        return (!((*this < second) || *this == second));
    }

    bool operator>=(const BigInteger& second)const {
        return (!(*this < second));
    }

    BigInteger &operator+=(const BigInteger& second) {
        if (sign == second.sign) {
            for (long long i = 0; i < second.size; i++) {
                if (size <= i) {
                    size++;
                    numbers.push_back(0);
                }
                numbers[i] += (second.numbers[i]);
            }
            for (long long i = 0; i < size - 1; i++) {
                numbers[i + 1] += numbers[i] / mod;
                numbers[i] %= mod;
            }
            while (numbers[size - 1] >= mod) {
                numbers.push_back(numbers[size - 1] / mod);
                numbers[size - 1] %= mod;
                size++;
            }
            return *this;
        }else{
            bool left_abs_more = ((sign && (-*this > second)) || (!sign && (*this > -second)));
            for (int i = 0; i < second.size; i++){
                if (size <= i) {
                    size++;
                    numbers.push_back(0);
                }
                if (left_abs_more) {
                    numbers[i] -= second.numbers[i];
                } else {
                    numbers[i] = -numbers[i] + second.numbers[i];
                }
            }

            if (!left_abs_more) {
                sign = !sign;
            }
            for (long long i = 0; i < size - 1; i++) {
                if (numbers[i] < 0) {
                    numbers[i] += mod;
                    numbers[i + 1]--;
                }
            }
            delete_zero();

            if (size == 1 && numbers[0] == 0) {
                sign = 0;
            }
            return *this;
        }
    }

    BigInteger &operator-=(const BigInteger& second) {
        *this += (-second);
        return *this;
    }

    BigInteger &operator*=(const BigInteger& second) {
        sign = (sign == second.sign ? 0 : 1);
        for (int i = 0; i < second.size + 1; i++) {
            numbers.push_back(0);
        }
        long long now_size = size;
        long long now_size_second = second.size;
        for (int i = now_size - 1; i >= 0; i--) {
            for (int j = now_size_second - 1; j >= 0; j--) {
                if (j == 0) {
                    numbers[i] *= second.numbers[j];
                    continue;
                }
                numbers[i + j] += numbers[i] * second.numbers[j];
            }
            for (int j = 0; j < now_size + now_size_second; j++) {
                numbers[j + 1] += (numbers[j] / mod);
                numbers[j] %= mod;
            }
        }

        size += second.size + 1;
        now_size += now_size_second + 1;
        for (int i = 0; i < now_size - 1; i++) {
            numbers[i + 1] += (numbers[i] / mod);
            numbers[i] %= mod;
        }
        delete_zero();
        if (size == 1 && numbers[0] == 0) {
            sign = 0;
        }
        return *this;
    }

    BigInteger &operator/=(const int del_old) {
        long long del = del_old;
        if (del < 0){
            sign = (-sign);
            del = -del;
        }
        for (long long i = size - 1; i >= 0; i--){
            if (i > 0){
                numbers[i - 1] += (numbers[i] % del) * mod;
            }
            numbers[i] /= del;
        }
        delete_zero();
        return *this;
    }

    BigInteger &operator/=(const BigInteger& divider) {
        bool new_sign = sign;
        if (divider.sign){
            new_sign = 1 - sign;
        }
        sign = 0;
        std::vector<long long> now;
        BigInteger check;
        check.numbers = divider.numbers;
        check.size = divider.size;
        check.sign = 0;
        for (long long i = size - 1; i >= (divider.size - 1); i--) {
            long long l = 0, r = 1e9;
            while (r - l > 1){
                long long mid = (l + r) / 2;
                if ((long long)(check.numbers.size()) != divider.size) {
                    check.numbers.pop_back();
                    check.size--;
                }
                for (int j = 0; j < divider.size; j++) {
                    check.numbers[j] = divider.numbers[j] * mid;
                }
                for (int j = 0; j < divider.size - 1; j++) {
                    check.numbers[j + 1] += (check.numbers[j] / mod);
                    check.numbers[j] %= mod;
                }
                if (check.numbers[divider.size - 1] >= mod) {
                    check.numbers.push_back(check.numbers[divider.size - 1] / mod);
                    check.numbers[divider.size - 1] %= mod;
                    check.size++;
                }
                if (comp_shift(check, i - divider.size + 1)) {
                    l = mid;
                }else{
                    r = mid;
                }
            }

            if (check.numbers.size() != divider.numbers.size()) {
                check.numbers.pop_back();
                check.size--;
            }
            for (long long j = 0; j < divider.size; j++) {
                check.numbers[j] = divider.numbers[j] * l;
            }
            for (long long j = 0; j < divider.size - 1; j++) {
                check.numbers[j + 1] += (check.numbers[j] / mod);
                check.numbers[j] %= mod;
            }
            if (check.numbers[divider.size - 1] >= mod) {
                check.numbers.push_back(check.numbers[divider.size - 1] / mod);
                check.numbers[divider.size - 1] %= mod;
                check.size++;
            }
            if (check.size > divider.size) {
                numbers[i + 1] -= check.numbers[divider.size];
            }
            for (long long j = i; j > i - divider.size; j--) {
                numbers[j] -= check.numbers[j - (i - divider.size) - 1];
            }
            for (long long j = 0; j < size - 1; j++) {
                while (numbers[j] < 0) {
                    numbers[j + 1]--;
                    numbers[j] += mod;
                }
                if (numbers[j] >= mod) {
                    numbers[j + 1] += (numbers[j] / mod);
                    numbers[j] %= mod;
                }
            }
            delete_zero();
            now.push_back(l);
        }
        if (now.size() == 0) {
            now.push_back(0);
        }
        for (int i = 0; i < (int)(now.size()) / 2; i++) {
            std::swap(now[i], now[int(now.size()) - i - 1]);
        }
        numbers = now;
        size = now.size();
        sign = new_sign;
        return *this;
    }
    BigInteger &operator%=(const BigInteger& divider) {
        BigInteger now = *this;
        now /= divider;
        now *= divider;
        *this -= now;
        return *this;
    }



    BigInteger operator++(int) {
        BigInteger copy(*this);
        if (!sign) {
            int i = 0;
            if (size == 0) {
                size++;
                numbers.push_back(0);
            }
            numbers[i]++;
            while (i < size - 1 && numbers[i] >= mod) {
                ++numbers[i + 1];
                numbers[i] -= mod;
            }
            if (numbers[size - 1] >= mod) {
                numbers.push_back(1);
                numbers[size - 1] -= mod;
            }
        }else {
            int i = 0;
            if (size == 1 && numbers[0] == 1) {
                numbers[0] = 0;
                sign = 0;
                return copy;
            }
            numbers[i]--;
            while (i < size - 1 && numbers[i] < 0) {
                --numbers[i + 1];
                numbers[i] += mod;
            }
        }

        return copy;
    }

    BigInteger& operator++() {
        if (!sign) {
            int i = 0;
            if (size == 0) {
                size++;
                numbers.push_back(0);
            }
            numbers[i]++;
            while (i < size - 1 && numbers[i] >= mod) {
                ++numbers[i + 1];
                numbers[i] -= mod;
            }
            if (numbers[size - 1] >= mod) {
                numbers.push_back(1);
                numbers[size - 1] -= mod;
            }
        }else {
            int i = 0;
            if (size == 1 && numbers[0] == 1) {
                numbers[0] = 0;
                sign = 0;
                return *this;
            }
            numbers[i]--;
            while (i < size - 1 && numbers[i] < 0) {
                --numbers[i + 1];
                numbers[i] += mod;
            }
        }

        return *this;
    }

    BigInteger operator--(int) {
        BigInteger copy(*this);
        if (sign) {
            int i = 0;
            numbers[i]++;
            while (i < size - 1 && numbers[i] >= mod) {
                ++numbers[i + 1];
                numbers[i] -= mod;
            }
            if (numbers[size - 1] >= mod) {
                numbers.push_back(1);
                numbers[size - 1] -= mod;
            }
        }else {
            int i = 0;
            if (size == 0) {
                size++;
                numbers.push_back(1);
                sign = 1;
                return copy;
            }
            if (size == 1 && numbers[0] == 0) {
                numbers[0] = 0;
                sign = -1;
                return copy;
            }
            numbers[i]--;
            while (i < size - 1 && numbers[i] < 0) {
                --numbers[i + 1];
                numbers[i] += mod;
            }
        }
        return copy;
    }

    BigInteger&operator--() {
        if (sign) {
            int i = 0;
            numbers[i]++;
            while (i < size - 1 && numbers[i] >= mod) {
                ++numbers[i + 1];
                numbers[i] -= mod;
            }
            if (numbers[size - 1] >= mod) {
                numbers.push_back(1);
                numbers[size - 1] -= mod;
            }
        }else {
            int i = 0;
            if (size == 0) {
                size++;
                numbers.push_back(1);
                sign = 1;
                return *this;
            }
            if (size == 1 && numbers[0] == 0) {
                numbers[0] = 0;
                sign = -1;
                return *this;
            }
            numbers[i]--;
            while (i < size - 1 && numbers[i] < 0) {
                --numbers[i + 1];
                numbers[i] += mod;
            }
        }

        return *this;
    }
    explicit operator bool() const {
        if (size >= 1 && (size > 1 || numbers[0] > 0)) {
            return 1;
        }
        return 0;
    }

    BigInteger operator - () const {
        BigInteger copy(*this);
        if (copy.size == 0 || (copy.size == 1 && numbers[0] == 0)) {
            return copy;
        }
        copy.sign = 1 - copy.sign;
        return copy;
    }
    std::string toString() const {
        std::string ans = "";
        long long now_size = size;
        if (size == 1 && numbers[0] == 0){
            ans += '0';
            return ans;
        }
        while (now_size > 0 && numbers[now_size - 1] == 0) {
            now_size--;
        }
        if (sign){
            ans += '-';
        }
        std::string cur_digit = "";
        long long div = 10;
        for (int i = 0; i < 9; i++) {
            cur_digit += char((numbers[now_size - 1] % div) / (div / 10) + '0');
            div *= 10;
        }
        for (int i = 8; i >= 0; i--) {
            if (cur_digit[i] =='0') {
                cur_digit.pop_back();
            }else{
                break;
            }
        }
        for (long long i = 0; i < (long long)(cur_digit.size()) / 2; i++) {
            std::swap(cur_digit[i], cur_digit[(long long)(cur_digit.size()) - 1 - i]);
        }
        ans += cur_digit;
        for (int i = now_size - 2; i >= 0; i--) {
            div = mod;
            cur_digit.clear();
            for (int j = 0; j < 9; j++){
                cur_digit += char(((numbers[i] % div) / (div / 10)) + '0');
                div /= 10;
            }
            ans += cur_digit;
        }
        return ans;
    }
    friend double to_double(const BigInteger& second);
    friend BigInteger gcd(const BigInteger &first, const BigInteger &second);
    friend std::istream& operator>>(std::istream& in, BigInteger& second);

};


BigInteger operator+(const BigInteger& first, const BigInteger& second) {
    BigInteger new_big;
    new_big = first;
    new_big += second;
    return new_big;
}

BigInteger operator-(const BigInteger& first,const BigInteger& second) {
    BigInteger new_big = first;
    new_big -= second;
    return new_big;
}


BigInteger operator*(const BigInteger& first,const BigInteger& second) {
    BigInteger new_big = first;
    new_big *= second;
    return new_big;
}

BigInteger operator/(const BigInteger& first,const BigInteger& second) {
    BigInteger new_big = first;
    new_big /= second;
    return new_big;
}

BigInteger operator%(const BigInteger& first,const BigInteger& second) {
    BigInteger new_big = first;
    new_big %= second;
    return new_big;
}

std::ostream& operator<<(std::ostream& out, const BigInteger& first) {
    std::string ans = first.toString();
    out << ans;
    return out;
}

std::istream& operator>>(std::istream& in, BigInteger& Big_in) {
    std::string str;
    in >> str;
    int id = 0;
    Big_in.sign = 0;
    if (str[id] == '-') {
        Big_in.sign = 1;
        id++;
    }
    Big_in.size = 0;
    Big_in.numbers.clear();
    long long col = 1;
    long long pos = 1;
    for (long long i = (long long)(str.size()) - 1; i >= id; i--) {
        if (Big_in.size < (col / 9) + std::min((long long)1, col % 9)) {
            Big_in.size++;
            pos = 1;
            Big_in.numbers.push_back(0);
        }
        Big_in.numbers[Big_in.size - 1] += pos * (long long)(str[i] - '0');
        col++;
        pos *= 10;
    }
    return in;
}


BigInteger gcd(const BigInteger &first, const BigInteger &second) {
    BigInteger now_first = first;
    BigInteger now_second = second;
    if (now_first < 0) {
        now_first = -now_first;
    }
    if (now_second < 0) {
        now_second = -now_second;
    }
    BigInteger ans = 1;
    while (now_second > 1 && now_first > 1 && (now_first != now_second)) {
        if (now_first.numbers[0] % 2 == 0 && now_second.numbers[0] % 2 == 0) {
            ans *= 2;
            now_first /= 2;
            now_second /= 2;
            continue;
        }
        if (now_first != 0 && now_first.numbers[0] % 2  == 0) {
            now_first /= 2;
            continue;
        }
        if (now_second != 0 && now_second.numbers[0] % 2 == 0) {
            now_second /= 2;
            continue;
        }
        if (now_first > now_second) {
            now_first -= now_second;
            now_first /= 2;
        } else {
            now_second -= now_first;
            now_second /= 2;
        }
    }
    if (now_first == 1 || now_second == 1) {
        return ans;
    }
    if (now_first > 0){
        return now_first * ans;
    }
    return now_second * ans;
}

double to_double(const BigInteger& now) {
    double ans = 0;
    double pow = 1;
    for (int i = 0; i < now.size; i++) {
        ans += pow * now.numbers[i];
        pow *= (1e9);
    }
    return ans;
}

class Rational {
private:
    BigInteger numerator;
    BigInteger denumerator;
    bool sign;

public:
    Rational() = default;
    Rational(const long long& number) {
        numerator = BigInteger(number);
        denumerator = BigInteger(1);
        if (number < 0) {
            sign = 1;
            numerator = (-numerator);
        }else {
            sign = 0;
        }
    }

    Rational(const BigInteger& number) {
        numerator = number;
        if (numerator < 0) {
            sign = 1;
            numerator = (-numerator);
        } else {
            sign = 0;
        }
        denumerator = BigInteger(1);
    }

    void swap(Rational& second) {
        std::swap(numerator, second.numerator);
        std::swap(denumerator, second.denumerator);
        std::swap(sign, second.sign);
    }

    Rational &operator=(const Rational& second) {
        if (&second == this) {
            return *this;
        }
        Rational copy(second);
        swap(copy);
        return *this;
    }

    bool operator==(const Rational& second) const {
        return (numerator == second.numerator && denumerator == second.denumerator && sign == second.sign);
    }

    bool operator!=(const Rational& second) const {
        return (!(*this == second));
    }

    bool operator<(const Rational& second) const {
        if (sign && second.sign == 0) {
            return true;
        }
        if (sign == 0 && second.sign) {
            return false;
        }
        BigInteger left = numerator * second.denumerator;
        BigInteger right = denumerator * second.numerator;
        if (sign == 0){
            return left < right;
        }else{
            return left > right;
        }
    }

    bool operator>(const Rational& second) const {
        return (!(*this < second || *this == second));
    }

    bool operator<=(const Rational& second) const {
        return (*this < second || *this == second);
    }

    bool operator>=(const Rational& second) const {
        return (!(*this < second));
    }

    Rational &operator+=(const Rational& second) {
        if (sign == second.sign) {
            numerator = numerator * second.denumerator + denumerator * second.numerator;
            denumerator *= second.denumerator;
            relax();
        }else{
            if (sign){
                if ((-*this) <= second) {
                    sign = 0;
                }
            }else{
                if (*this <= (-second)) {
                    sign = 1;
                }
            }
            numerator = denumerator * second.numerator - numerator * second.denumerator;
            denumerator *= second.denumerator;
            if (numerator < 0) {
                numerator = -numerator;
            }
            relax();
        }
        if (numerator == 0) {
            denumerator = 1;
            sign = 0;
        }
        return *this;
    }

    Rational &operator-=(const Rational& second) {
        *this += (-second);
        return *this;
    }

    Rational &operator*=(const Rational& second) {
        if (second.sign) {
            sign = 1 - sign;
        }
        numerator *= second.numerator;
        denumerator *= second.denumerator;
        relax();
        return *this;
    }

    Rational &operator/=(const Rational& second) {
        if (second.sign) {
            sign = 1 - sign;
        }
        numerator *= second.denumerator;
        denumerator *= second.numerator;
        relax();
        return *this;
    }

    Rational operator - () const {
        Rational copy(*this);
        copy.sign = 1 - copy.sign;
        return copy;
    }

    void relax() {
        if (numerator == 0){
            denumerator = 1;
            sign = 0;
            return;
        }
        BigInteger g = gcd(numerator, denumerator);
        numerator /= g;
        denumerator /= g;
    }

    std::string toString() const {
        std::string answer = "";
        if (sign){
            answer += '-';
        }
        answer += numerator.toString();
        if (denumerator != 1){
            answer += '/';
            answer += denumerator.toString();
        }
        return answer;
    }

    std::string asDecimal(size_t size) {
        BigInteger now = numerator;
        BigInteger ans = now / denumerator;
        std::string answer = "";
        if (sign){
            answer += '-';
        }
        answer += ans.toString();
        if (size > 0) {
            answer += '.';
            now %= denumerator;
            std::string now_s;
            for (size_t i = 0; i < size; i++) {
                now *= 10;
                ans = now / denumerator;
                now_s = ans.toString();
                now %= denumerator;
                answer += now_s[int(now_s.size()) - 1];
            }

        }
        return answer;
    }

    explicit operator double() const {
        double ans = 0;
        BigInteger now = numerator;
        BigInteger tmp;
        double pow = 1;
        for (int i = 0; i < 10; i++) {
            tmp = now;
            now /= denumerator;
            ans += to_double(now) * pow;
            now = tmp - now * denumerator;
            now *= 1e9;
            pow /= 1e9;
        }
        if (sign) {
            ans *= -1;
        }
        return ans;
    }
    friend std::istream& operator>>(std::istream& in, Rational& Big_input);
};

std::ostream& operator<<(std::ostream& out, const Rational& now){
    std::string string_now = now.toString();
    out << string_now;
    return out;
}

std::istream& operator>>(std::istream& in, Rational& ans) {
    BigInteger now;
    in >> now;
    ans.numerator = now;
    ans.denumerator = 1;
    ans.sign = 0;
    if (ans.numerator < 0){
        ans.sign = 1;
        ans.numerator = (-ans.numerator);
    }
    return in;
}

Rational operator+(const Rational& first, const Rational& second){
    Rational new_big;
    new_big = first;
    new_big += second;
    return new_big;
}

Rational operator-(const Rational& first,const Rational& second){
    Rational new_big = first;
    new_big -= second;
    return new_big;
}


Rational operator*(const Rational& first,const Rational& second){
    Rational new_big = first;
    new_big *= second;
    return new_big;
}

Rational operator/(const Rational& first,const Rational& second){
    Rational new_big = first;
    new_big /= second;
    return new_big;
}



