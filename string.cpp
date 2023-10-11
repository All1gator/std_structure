#include <iostream>
#include <cstring>

class String{
private:
    size_t size = 0;
    size_t buffer_size = 0;
    char* buffer = nullptr;


public:
    String() = default;

    String(size_t n, char c): size(n), buffer_size(2 * n), buffer(new char[buffer_size]) {
        memset(buffer, c, n);
    }

    String(const String& str): String(str.size, '\0') {
        memcpy(buffer, str.buffer, size);
    }

    String(const char c[]): size(strlen(c)), buffer_size(strlen(c) * 2), buffer(new char[buffer_size]) {
        memcpy(buffer, c, size);
    }

    void swap(String& str) {
        std::swap(buffer, str.buffer);
        std::swap(size, str.size);
        std::swap(buffer_size, str.buffer_size);
    }

    String &operator=(const String& str) {
        if (&str == this){
            return *this;
        }
        String copy(str);
        swap(copy);
        return *this;
    }

    bool operator==(const String& str) const {
        if (size != str.size) {
            return 0;
        }
        for (size_t i = 0; i < size; i++) {
            if (buffer[i] != str.buffer[i]) {
                return 0;
            }
        }
        return 1;
    }

    bool operator!=(const String& str) const {
        if (*this == str) {
            return 0;
        }else {
            return 1;
        }
    }

    String& operator+=(const String& str) {
        size_t old_size = str.size;
        for (size_t i = 0; i < old_size; i++) {
            push_back(str[i]);
        }
        return *this;
    }

    String& operator+=(const char& symbol) {
        push_back(symbol);
        return *this;
    }

    size_t find(const String& str) const {
        size_t index = size;
        if (size < str.size) {
            return index;
        }
        for (size_t i = 0; i <= size - str.size; ++i) {
            bool is_find = 1;
            for (size_t j = i; j < i + str.size; ++j) {
                if (buffer[j] != str.buffer[j - i]) {
                    is_find = 0;
                }
            }
            if (is_find){
                index = i;
                break;
            }
        }
        return index;
    }

    size_t rfind(const String& str) const {
        size_t index = size;
        if (size < str.size) {
            return index;
        }
        for (int i = int(size) - int(str.size); i >= 0; --i) {
            bool is_find = 1;
            for (size_t j = i; j < i + str.size; ++j) {
                if (buffer[j] != str.buffer[j - i]) {
                    is_find = 0;
                }
            }
            if (is_find) {
                index = i;
                break;
            }
        }
        return index;
    }

    String substr(int start, int count) const {
        String new_String = String(count, '\0');
        memcpy(new_String.buffer, (buffer + start), count);
        return new_String;
    }

    bool empty() const {
        return size == 0;
    }

    void clear() {
        size = 0;
        buffer_size = 0;
        delete[] buffer;
        buffer = nullptr;
    }

    ~String() {
        delete[] buffer;
    }

    size_t length() const {
        return size;
    }

    char& operator[](int index) {
        return buffer[index];
    }

    char operator[](int index) const {
        return buffer[index];
    }

    char& front() {
        return buffer[0];
    }

    char front() const {
        return buffer[0];
    }

    char &back() {
        return buffer[size - 1];
    }

    char back() const {
        return buffer[size - 1];
    }

    void push_back(char new_symbol) {
        if (buffer_size > size) {
            buffer[size] = new_symbol;
            ++size;
        } else {
            buffer_size *= 2;
            ++buffer_size;
            char* new_buffer = new char[buffer_size];
            memcpy(new_buffer, buffer, size);
            buffer = new_buffer;
            buffer[size] = new_symbol;
            ++size;

        }
    }

    void pop_back() {
        --size;
        if (buffer_size >= 4 * size ) {
            char* new_buffer = new char[buffer_size / 2];
            memcpy(new_buffer, buffer, size);
            buffer = new_buffer;
            buffer_size /= 2;
        }

    }
};


String operator+(const String& first, const String& second) {
    String new_String = String(first);
    new_String += second;
    return new_String;
}

String operator+(const char& first, const String& second) {
    String new_String;
    new_String.push_back(first);
    new_String += second;
    return new_String;
}

String operator+(const String& first, const char& second) {
    String new_String = String(first);
    new_String.push_back(second);
    return new_String;
}

std::ostream& operator<<(std::ostream& out, const String& str) {
    for (size_t i = 0; i < str.length(); i++) {
        out << str[i];
    }
    return out;
}

std::istream& operator>>(std::istream& in, String& str) {
    str.clear();
    if (!in.eof()) {
        char symbol;
        symbol = in.get();
        while ((isspace(symbol) || symbol == '\n') && (!in.eof())) {
            symbol = in.get();
        }
        while ((!isspace(symbol) && symbol != '\n') && (!in.eof())) {
            str.push_back(symbol);
            symbol = in.get();
        }
    }
    return in;
}
