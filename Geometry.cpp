#include <iostream>
#include <cstring>
#include <vector>
#include <string>
#include <cassert>
#include <cmath>

double eps = (1.0) / 1e4;
double PI = 3.1415926535;

bool double_equal(double first, double second) {
    return (abs(first - second) < eps ? 1 : 0);
}


struct Point {
    double x, y;
    Point() = default;
    Point(double now_x, double now_y): x(now_x), y(now_y) {};

    bool operator==(const Point& second) const {
        return (double_equal(x, second.x)&& double_equal(y, second.y));
    }

    bool operator!=(const Point& second) const {
        return (!(*this == second));
    }

    Point &operator+=(const Point& second) {
        x += second.x;
        y += second.y;
        return *this;
    }



    Point &operator-=(const Point& second) {
        x -= second.x;
        y -= second.y;
        return *this;
    }

    Point &operator*=(const double& second) {
        x *= second;
        y *= second;
        return *this;
    }

    Point &operator/=(const double& second) {
        x /= second;
        y /= second;
        return *this;
    }

};

Point operator+(const Point& first, const Point& second) {
    Point now = first;
    now += second;
    return now;
}

Point operator-(const Point& first, const Point& second) {
    Point now = first;
    now -= second;
    return now;
}

Point operator*(const Point& first, const double& second) {
    Point now = first;
    now *= second;
    return now;
}

Point operator/(const Point& first, const double& second) {
    Point now = first;
    now /= second;
    return now;
}

double vector_mul(const Point& first, const Point& second) {
    return (first.x) * (second.y) - (first.y) * (second.x);
}

class Line{
private:
    Point r_zero;
    Point r;
public:

    Line(const Point& first, const Point& second) {
        r_zero = first;
        r = second - first;
    }

    Line(const double& k, const double& b) {
        r.x = 1;
        r.y = k;
        r_zero.x = 0;
        r_zero.y = b;
    }

    Line(const Point& first, const double& k) {
        r.x = 1;
        r.y = k;
        r_zero = first;
    }

    Line normal(const Point& first) const {
        Line now(first, first);
        Point revers = r;
        std::swap(revers.x, revers.y);
        revers.x = -revers.x;
        now.update(revers, first);
        return now;
    }


    bool operator==(const Line& second) const {
        return (double_equal(vector_mul(r, second.r), 0) && double_equal(vector_mul(r_zero - second.r_zero, r), 0));
    }

    bool operator!=(const Line& second) const {
        return (!(*this == second));
    }

    void update( Point& first,const Point& second) {
        r = first;
        r_zero = second;
    }

    Point intersection(const Line& second) const {
        Point ans = r * (vector_mul(second.r_zero, second.r));
        ans -= second.r * (vector_mul(r_zero, r));
        ans /= vector_mul(r, second.r);
        return ans;
    }

    bool incle(const Point& a) const {
        return (double_equal(vector_mul(r, a - r_zero), 0));
    }



};

bool sign_vector_mul(const Point& first, const Point& second, const Point& third) {
    return (((second.x - first.x) * (third.y - second.y) - (second.y - first.y) * (third.x - second.x)) > 0 ? 1 : 0);
}


double distant (const Point& first, const Point& second) {
    return sqrt((first.x - second.x) * (first.x - second.x) + (first.y - second.y) * (first.y - second.y));
}


class Shape {
protected:
    Point focus_1, focus_2;
    double a;
    double extr;
    double b;
    std::vector<Point> verticles;

public:
    Shape() = default;
    Shape(const Shape& now) {
        focus_1 = now.focus_1;
        focus_2 = now.focus_2;
        a = now.a;
        extr = now.extr;
        b = now.b;
        verticles = now.verticles;
        relax();
    }

    Shape &operator=(const Shape& now) {
        focus_1 = now.focus_1;
        focus_2 = now.focus_2;
        a = now.a;
        extr = now.extr;
        b = now.b;
        verticles = now.verticles;
        return *this;

    }

    virtual ~Shape() = default;
    virtual double area() const = 0;
    virtual double perimeter() const = 0;

    bool isSimilarTo(const Shape& another) const{
        if (verticles.size() == 0 && another.verticles.size() == 0) {
            return (a / another.a == b / another.b);
        }
        if (verticles.size() != another.verticles.size()) {
            return false;
        }
        bool yes = 0;
        int n = verticles.size();

        for (int i = 0; i < n; i++) {
            double k = distant(verticles[1], verticles[0]);
            k /= distant(another.verticles[i], another.verticles[(i + 1) % n]);
            yes = 1;
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n; k++) {
                    if (!double_equal(distant(verticles[j], verticles[(j + 1) % n]),
                                      distant(verticles[(j + i) % n], verticles[(j + 1 + i) % n]))) {
                        yes = 0;
                    }
                }
            }
            if (yes) {
                return 1;
            }
        }

        return false;
    }

    bool  isCongruentTo(const Shape& any) const {
        if (verticles.size() == 0 && any.verticles.size() == 0) {
            return (double_equal(a,any.a) && double_equal(b, any.b));
        }

        if (verticles.size() != any.verticles.size()) {
            return false;
        }

        int n = int(verticles.size());
        /*std::cerr << "ALARMMMMMMMM" << '\n';
        std::cerr << verticles.size() <<  '\n';
        for (int i = 0; i < int(verticles.size()); i++) {
            std::cerr << distant(verticles[i], verticles[prev(i)]) << '\n';
        }
        std::cerr << '\n';

        std::cerr << verticles.size() <<  '\n';
        for (int i = 0; i < int(verticles.size()); i++) {
            std::cerr << distant(any.verticles[i], any.verticles[prev(i)]) << '\n';
        }
        std::cerr << '\n';
        std::cerr << '\n';

        std::cerr << "FIRST" << '\n';
        std::cerr << verticles.size() << '\n';


        for (int i = 0; i < n; i++) {
            std::cerr << verticles[i].x << ' ' << verticles[i].y << '\n';
        }
        std::cerr << '\n';
        for (int i = 0; i < n; i++) {
            std::cerr << any.verticles[i].x << ' ' << any.verticles[i].y << '\n';
        }
        std::cerr << '\n';*/

        for (int i = 0; i < n; i++) {
            Point now = any.verticles[next(i)] - any.verticles[i];
            std::cerr << now.x << ' ' << now.y  << ' ' << distant(any.verticles[next(i)], any.verticles[i]) << '\n';
        }
        std::cerr << '\n';

        for (int i = 0; i < n; i++){
            bool is_equ = 1;
            for (int j = 0; j < n; j++) {
                Point first = verticles[j] - verticles[prev(j)];
                Point second = verticles[next(j)] - verticles[j];
                double dist_left = distant(verticles[next(j)], verticles[j]);
                double left = vector_mul(first, second);

                first = any.verticles[(j + i) % n] - any.verticles[prev((j + i) % n)];
                Point second_2 = any.verticles[next((j + i )%n)] - any.verticles[(j + i) % n];
                double dist_right = distant(any.verticles[next((j + i)%n)], any.verticles[(j + i) % n]);
                double right = vector_mul(first, second_2);
                if (!double_equal(left, right) || !double_equal(dist_left, dist_right)) {
                    is_equ = 0;
                    //std::cerr << i << ' ' << j << ' ' << left << ' ' << right << '\n';
                    break;
                }
            }
            if (is_equ) {
                std::cerr << 1 << '\n';
                return true;
            }
            is_equ = 1;
            for (int j = 0; j < n; j++) {
                Point first = verticles[j] - verticles[prev(j)];
                Point second = verticles[next(j)] - verticles[j];
                double dist_left = distant(verticles[next(j)], verticles[j]);
                double left = vector_mul(first, second);

                first = any.verticles[(-j + i + n) % n] - any.verticles[next((-j + i + n) % n)];
                Point second_2 = any.verticles[prev((-j + i + n)%n)] - any.verticles[(-j + i + n) % n];
                double dist_right = distant(any.verticles[prev((-j + i + n)%n)], any.verticles[(-j + i + n) % n]);
                double right = vector_mul(second_2, first);
                if (!double_equal(left, right) || !double_equal(dist_left, dist_right)) {
                    is_equ = 0;
                    //std::cerr << i << ' ' << j << ' ' << left << ' ' << right << '\n';
                    break;
                }
            }

            if (is_equ) {
                std::cerr << 1 << '\n';
                return true;
            }
        }
        std::cerr << 0 << '\n';
        return false;
    }

    bool containsPoint(const Point& point) const{
        if (verticles.size() == 0) {
            return (distant(point, focus_1) + distant(point, focus_2) <= 2 * a );
        }
        Line beam(point, point);
        Point now;
        now.x = 1337;
        now.y = 79;
        beam.update(now, point);
        int cnt = 0;
        int n = verticles.size();
        for (int i = 0; i < n; i++) {
            if (beam.incle(verticles[i])) {
                cnt++;
            }
        }
        int col = 0;
        for (int i = 0; i < n; i++) {
            Line now( verticles[(i + 1) % n], verticles[i]);
            Point first = now.intersection(beam);
            Point now_vec = first - point;
            if (now_vec.x  < 0) {
                continue;
            }
            if (((verticles[(i + 1) % n].x <= point.x && point.x <= verticles[(i) % n].x) ||
                (verticles[(i) % n].x <= point.x && point.x <= verticles[(i + 1) % n].x)) &&
                ((verticles[(i + 1) % n].y <= point.y && point.y <= verticles[(i) % n].y) ||
            (verticles[(i) % n].y <= point.y && point.y <= verticles[(i + 1) % n].y))) {
                col++;
            }
        }
        col -= cnt;
        return (col % 2);

    }



    Shape& rotate(const Point& center, double angle) {


        Point now = focus_1 - center;
        Point new_focus_1;
        new_focus_1.x = now.x * cos(angle * PI / 180) - now.y * sin(angle * PI / 180);
        new_focus_1.y = now.x * sin(angle * PI / 180) + now.y * cos(angle * PI / 180);
        focus_1 = new_focus_1 + center;

        now = focus_2 - center;
        Point new_focus_2;
        new_focus_2.x = now.x * cos(angle * PI / 180) - now.y * sin(angle * PI / 180);
        new_focus_2.y = now.x * sin(angle * PI / 180) + now.y * cos(angle * PI / 180);
        focus_2 = new_focus_2 + center;

        int n = verticles.size();
        for (int i = 0; i < n; i++) {
            now = verticles[i] - center;
            new_focus_1.x = now.x * cos(angle * PI / 180) - now.y * sin(angle * PI / 180);
            new_focus_1.y = now.x * sin(angle * PI / 180) + now.y * cos(angle * PI / 180);
            verticles[i] = new_focus_1 + center;
        }

        relax();

        return *this;

    }

    size_t next(int i) const{
        return ((i + 1) % (verticles.size()));
    }
    size_t prev(int i) const{
        return ((i - 1 >= 0) ? (i - 1) : (verticles.size() - 1));
    }

    Shape& reflect(const Line& axis) {
        //std:: cerr << "Reflect Line" << ' ' << number << ' ' << '\n';
        if (verticles.size() == 0) {
            Line now = axis.normal(focus_1);
            Point my = now.intersection(axis);
            Point add = my - focus_1;
            add *= 2;
            focus_1 += add;
            now = axis.normal(focus_2);
            my = now.intersection(axis);
            add = my - focus_2;
            add *= 2;
            focus_2 += add;
        }else{
            int n = verticles.size();
            Point my;
            Point add;
            Line now(my, my);
            for (int i = 0; i < n; i++) {
                now = axis.normal(verticles[i]);
                my = now.intersection(axis);
                add = my - verticles[i];
                add *= 2;
                verticles[i] += add;
            }
        }
        relax();
        return *this;
    }

    Shape& scale(const Point& center, const double& k) {
        //std:: cerr << "scale" << ' ' << number << ' ' << center.x << ' ' << center.y << ' ' << k << '\n';
        if (verticles.size() == 0) {
            Point now = focus_1 - center;
            now *= k;
            focus_1 = now + center;
            now = focus_2 - center;
            now *= k;
            focus_2 = now + center;
            a *= k;
            b *= k;
        }else{
            int n = verticles.size();
            Point now;
            for (int i = 0; i < n; i++){
                now = verticles[i] - center;
                now *= k;
                verticles[i] = now + center;
            }
        }

        relax();
        return *this;
    }

    Shape& relax() {
        if (verticles.size() == 0) {
            return *this;
        }
        int n = verticles.size();
        double minn = verticles[0].x;
        int start = 0;
        for (int i = 1; i < n; i++) {
            if (verticles[i].x < minn) {
                minn = verticles[i].x;
                start = i;
            }
        }
        std::vector<Point> copy = verticles;
        if (sign_vector_mul(verticles[prev(start)], verticles[start], verticles[next(start)])) {
            for (int i = 0; i <n; i++) {
                copy[i] = verticles[(i + start) % n];
            }
        }else{
            for (int i = 0; i < n; i++) {
                copy[i] = verticles[(start - i + n) % n];
            }
        }
        verticles = copy;
        return *this;
    }

    Shape& reflect(const Point& center) {
        //std:: cerr << "Reflectpoint" << ' ' << number << ' ' << center.x << ' ' << center.y  << '\n';
        int n = verticles.size();
        if (n > 0) {
            Point now;
            Point new_focus_1;
            for (int i = 0; i < n; i++) {
                now =  center - verticles[i];
                now *= 2;
                verticles[i] += now;
            }
            relax();
        }else {
            Point now;
            now =  center - focus_1;
            now *= 2;
            focus_1 += now;
            now =  center - focus_2;
            now *= 2;
            focus_2 += now;
        }
        relax();
        /*
        std::cerr << verticles.size() <<  '\n';
        for (int i = 0; i < int(verticles.size()); i++) {
            std::cerr << verticles[i].x << ' ' << verticles[i].y << '\n';
        }
        std::cerr << '\n';*/
        return *this;
    }

    bool operator==(const Shape& first) const {
        if (verticles.size() == 0 && first.verticles.size() == 0){
            return ((focus_1 == first.focus_1 && focus_2 == first.focus_2 && double_equal(a, first.a))
                    ||((focus_2 == first.focus_1 && focus_1 == first.focus_2 && double_equal(a, first.a))));
        }else  {

            if (first.verticles.size() != verticles.size()) {
                return 0;
            }
            int size = verticles.size();
            int start = -1;
            for (int i = 0; i < size; i++){
                if (first.verticles[i] == verticles[0]) {
                    start = i;
                    break;
                }
            }
            if (start == -1) {
                return 0;
            }
            for (int i = 0; i < size; i++) {
                if (first.verticles[(i + start) % size] != verticles[i]) {
                    return 0;
                }
            }
            return 1;
        }
    }
    bool operator!=(const Shape& first) const {
        return !(first == *this);
    }

};



class Ellipse : public Shape {
protected:


public:
    Ellipse() = default;

    Ellipse(const Point& first, const Point& second, const double& distance) {
        focus_1 = first;
        focus_2 = second;
        a = distance / 2;
        extr = (distant(first, second) / 2) / a;
        b = a * sqrt(1 - extr * extr);

    }





    std::pair<Point,Point> focuses() const{
        return std::make_pair(focus_1, focus_2);
    }

    double eccentricity() const {
        return extr;
    }



    double area() const override {
        return a * b * PI;
    }

    double perimeter() const override {
        return PI * (3 * (a + b) - sqrt((3 * a + b) * (3 * b + a)));
    }

};

class Circle  : public Ellipse {
protected:

public:
    Circle() = default;

    Circle(const Point& center, const double& distance) {
        focus_1 = center;
        focus_2 = center;
        a = distance;
        extr = 0;
        b = a;
    }
    Point center() const{
        return focus_1;
    }


    double radius() const {
        return a;
    }

};



class Polygon : public Shape{
protected:


public:

    Polygon() = default;

    Polygon(std::vector<Point>& ver) {
        verticles = ver;
        relax();
    }


    double perimeter() const override{
        if (verticles.size() <= 1) {
            return 0;
        }
        double ans = distant(verticles[0], verticles[verticles.size() - 1]);
        for (size_t i = 0; i < verticles.size() - 1; i++) {
            ans += distant(verticles[i], verticles[i + 1]);
        }
        return ans;
    }


    virtual double area() const override {
        if (verticles.size() < 3) {
            return 0;
        }
        double ans = 0;
        Point left, right;
        for (size_t i = 1; i < verticles.size() - 1; i++) {
            left = verticles[i] - verticles[0];
            right = verticles[i + 1] - verticles[0];
            ans += (vector_mul(left, right));
        }
        ans /= 2;
        return ans;
    }



    size_t verticesCount() const {
        return verticles.size();
    }

    std::vector<Point> getVertices() const {
        return verticles;
    }

    bool isConvex() const{
        if (verticles.size() < 3) {
            return true;
        }
        bool ans = sign_vector_mul(verticles[prev(0)], verticles[0], verticles[next(0)]);
        for (int i = 1; i < int(verticles.size()); i++) {
            if (ans != sign_vector_mul(verticles[prev(i)], verticles[i], verticles[next(i)])) {
                return false;
            }
        }
        return true;
    }

    void unknown_members(std::vector<Point>& points)const;

    template<typename T, typename... Args>
    void unknown_members(std::vector<Point>& points, T start, Args... args)const;

    template<typename... Args>
    Polygon(Args... args);
};

void Polygon::unknown_members(std::vector<Point>& points) const {
    if (points.size()) {
        return;
    }
}

template<typename T, typename... Args>
void Polygon::unknown_members(std::vector<Point>& points, T start, Args... args) const {
    points.push_back(Point(start));
    unknown_members(points, args...);
}

template<typename... Args>
Polygon::Polygon(Args... args) {
    std::vector<Point> points;
    unknown_members(points, args...);
    verticles = points;
    relax();
    points.clear();
}




class Rectangle : public Polygon {
protected:


public:
    Rectangle() = default;

    Rectangle(const Point& first, const Point& second, const double k) {
        verticles.push_back(first);
        Point vec = second;
        vec -= first;
        double dist = distant(first, second);
        double t = dist / (1.0/k + k);
        Point now = first;
        Point vec_copy = vec;
        vec_copy /= (t/k + t * k);

        vec = vec_copy;
        vec *= (t / k);
        now += vec;
        vec_copy *= t;
        std::swap(vec_copy.x, vec_copy.y);
        vec_copy.x = -vec_copy.x;
        now += vec_copy;

        Point new_vec = first;
        new_vec -= now;
        new_vec += second;
        verticles.push_back(new_vec);
        verticles.push_back(second);
        verticles.push_back(now);
        relax();
    }

    virtual double area() const override {
        double first = distant(verticles[1], verticles[0]);
        double second = distant(verticles[3], verticles[0]);
        return first * second;
    }

    Point center() const{
        Point now = verticles[0] + verticles[2];
        now /= 2;
        return (now);
    }

    std::pair<Line, Line> diagonals() const{
        return {(Line(verticles[0], verticles[2])), Line(verticles[1], verticles[3])};
    }

};

class Square : public Rectangle {
protected:


public:
    Square() = default;

    Square(const Point& first, const Point& second) {
        verticles.push_back(first);
        Point now;
        now.x = second.x;
        now.y = first.y;
        verticles.push_back(now);
        verticles.push_back(second);
        now.x = first.x;
        now.y = second.y;
        verticles.push_back(now);
        relax();
    }

    double area() const override {
        double first = distant(verticles[1], verticles[0]);
        return first * first;
    }

    Circle inscribedCircle() const {
        Circle now(center(), distant(verticles[0], verticles[1]));
        return now;
    }

    Circle circumscribedCircle() const {
         Circle now(center(), distant(verticles[0], verticles[2]));
         return now;
    }

};

class Triangle : public Polygon {
protected:

public:
    using Polygon::Polygon;
    double area() const override {
        Point first = verticles[2] - verticles[0];
        Point second = verticles[1] - verticles[0];
        return abs(vector_mul(first, second)) / 2;
    }

    Circle inscribedCircle() const {
        double k = distant(verticles[0], verticles[1]) / distant(verticles[0], verticles[2]);
        double y = distant(verticles[1], verticles[2]) / (1 + k);
        Point vec_cb = verticles[1] - verticles[2];
        vec_cb /= distant(verticles[1], verticles[2]);
        vec_cb *= y;
        vec_cb += verticles[2];
        Line first(verticles[0], vec_cb);

        k = distant(verticles[2], verticles[1]) / distant(verticles[0], verticles[2]);
        double x = distant(verticles[0], verticles[1]) / (1 + k);
        Point vec_ab = verticles[1] - verticles[0];
        vec_ab /= distant(verticles[0], verticles[1]);
        vec_ab *= x;
        vec_ab += verticles[0];
        Line second(verticles[2], vec_ab);

        Point center = first.intersection(second);
        double r = 2 * area() / (distant(verticles[0], verticles[1]) + distant(verticles[0], verticles[2]) + distant(verticles[1], verticles[2]));
        Circle ans(center, r);
        return ans;
    }

    Circle circumscribedCircle() const {
        Point x = (verticles[0] + verticles[1]) / 2;
        Point now = verticles[1] - verticles[0];
        std::swap(now.x, now.y);
        now.x = -now.x;
        Point y = x + now;
        Line first(x, y);

        x = (verticles[0] + verticles[2]) / 2;
        now = verticles[2] - verticles[0];
        std::swap(now.x, now.y);
        now.x = -now.x;
        y = x + now;
        Line second(x, y);
        Point center = first.intersection(second);

        double r = distant(verticles[0], verticles[1]) * distant(verticles[0], verticles[2]) * distant(verticles[1], verticles[2]);
        r /= (4 *area());
        return Circle(center, r);
    }

    Point centroid() {
        return (verticles[0] + verticles[1] + verticles[2]) / 3;
    }

    Point orthocenter() {
        Line ac(verticles[0], verticles[1]);
        Line first = ac.normal(verticles[2]);
        Line ab(verticles[0], verticles[2]);
        Line second = ab.normal(verticles[1]);
        Point ans = first.intersection(second);
        std::cerr << ans.x << ' ' << ans.y << '\n';
        return first.intersection(second);
    }
    Line EulerLine() {
        return Line(centroid(), orthocenter());
    }

    Circle ninePointsCircle() {
        Circle cir = circumscribedCircle();
        cir.scale(orthocenter(), 0.5);
        return cir;
    }

};

