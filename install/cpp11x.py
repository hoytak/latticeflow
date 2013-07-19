from configure import chooseOption, filterOptions, addOption

using_cxx11 = False
cxx11_flag = ""

def setupCXX11(ctx):

    global using_cxx11, cxx11_flag

    using_cxx11 = True
    cxx11_flag = chooseOption(ctx, '-std=c++11', '-std=c++0x')


def checkCXX11features(ctx, required_features):

    def check(name, fragment, execute = False):

        ctx.check_cxx(
            msg = ("C++11x: Checking %s support for %s"
                   % ("required" if name in required_features else "optional", name)),
            define_name = 'CXX_SUPPORT_' + name.replace(' ', '_').upper(),
            fragment    = fragment,
            execute     = execute,
            mandatory   = (name in required_features) )

    required_features = set(required_features)

    ctx.env.stash()
    addOption(ctx, cxx11_flag)

    ############################################################

    try:

        check('auto',
              'int main() { auto x = 0; return x; }')

        check('lamda',
              """
              #include <algorithm>
              int main() {
                double x[2] = {1,2};
                std::sort(x, x + 2, [](double a, double b){return b < a;});
                return 0;
              }
              """)

        check('decltype',
              'int main() { int x = 0; decltype(x) y = x; return y; }')

        check('nullptr',
              """
              int main() { int *x = nullptr; int y = 0;
              x = &y; return *x;}
              """)

        check('rvalue reference',
              """
              int& rvalueTest(int&& x) { return x; }
              int main() { return 0; }
              """
              )

        check("static_assert",
              """
              int main() {
              static_assert(sizeof(int) == sizeof(int), \"test\");
              return 0;
              }
              """
              )

        check('strongly typed enums',
              """
              enum class TestEnum : unsigned char {kTestEnum1, kTestEnum2};
              int main() { return 0; }
              """
              )

        check('exception_ptr',
              """
              #include <stdexcept>

              std::exception_ptr epTest(std::exception_ptr eptr) {
                return eptr;
              }

              int main() { return 0; }
              """
              )

        check('__func__',
              """
              int main() {
                const char * test = __func__;
                test = 0;
                return (test != 0);
              }
              """
              )

        check('__FILE__',
              """
              int main() {
                char test[] = __FILE__;
                test[0] = 0;
                return test[0];
              }
              """
              )

        check('__LINE__',
              """
              int main() {
                int test = __LINE__;
                test = 0; return (test != 0);
              }
              """
              )

        check('_Pragma',
              """
              int main() {
                _Pragma ( \"pack()\" )
                return 0; }
              """
              )

        check('extended sizeof',
              """
              struct StructTest { int member; };
              int main() { return sizeof(StructTest::member); }
              """
              )

        check('trailing return types',
              """
              auto f() -> int;
              int main() { return 0; }
              """
              )

        check('right angle brackets',
              """
              #include <vector>
              typedef std::vector<std::vector<int>> Test;
              int main() { return 0; }
              """
              )

        check('extended friend declarations',
              """
              class Test { friend class TestTwo; };

            int main() { return 0; }
            """
            )

        check('range based for loop',
              """
            int main() {
                int array[6] = { 0, 1, 2, 3, 4, 5 };
                for (int& x : array) x *= 2;
                return array[0];
                }
            """
            )

        check('non-static initializers',
              """
            class Test { public: int a = 7; };
            int main() { return 0; }
            """
            )

        check('char16 type',
              "int main() { char16_t t = 0; return t; }"
              )

        check('char32 type',
              "int main() { char32_t t = 0; return t; }"
              )

        check('long long',
              'int main() { long long test = 0; return test;}'
              )

        check('initializer list',
              """
            #include <initializer_list>
            #include <vector>

            struct Test {
                std::vector<int> v;
                Test(std::initializer_list<int>l) : v(l){}; };

            int main() {
                Test *t = new Test({0, 1, 2});
                return !t;
            }
            """
            )

        check('variadic templates',
              """
            template<typename ...Args> struct Test {
                static const int size = sizeof...(Args);};
            int main() {
                Test<int, int>* t = new Test<int, int>;
                return !t; }
            """)


        check('SFINAE',
              """
            template <int I> struct A {};
            char xxx(int);char xxx(float);

            template <class T> A<sizeof(xxx((T)0))> f(T){A<1> b; return b;}

            int main() {f(4); return 0;}
            """
            )

        check('variadic macro support',
              """
            #define test(...) __VA_ARGS__
            int main() { test(return 0); return 0; }
            """,
            execute = True
          )

        check('constexpr',
              """
            struct A {
                static constexpr int x = 0;
                };
            int main() { return A::x; }
            """)

        check('template arguments',
              """
              #include <string>
              #include <map>
              #include <unordered_map>

              template <template <class...> class T>
              void test()
              {
                 T<std::string, int> x;
              }

              int main()
              {
                 test<std::map>();
                 test<std::unordered_map>();
              }
              """)


        check('template aliasing',
              """
              #include <map>
              template <class A>
              using str_map = map<A, str>

              int main()
              {
                  str_map<int> s;
              }
              """)
    ########################################
    # Tear down

    finally:

        ctx.env.revert()
