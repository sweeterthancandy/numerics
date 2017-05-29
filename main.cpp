#include <iostream>
#include <vector>
#include <functional>
#include <cassert>

#include <boost/lexical_cast.hpp>
#include <boost/preprocessor.hpp>

#define PRINT_SEQ_detail(r, d, i, e) do{ std::cout << ( i ? ", " : "" ) << BOOST_PP_STRINGIZE(e) << " = " << (e); }while(0);
#define PRINT_SEQ(SEQ) do{ BOOST_PP_SEQ_FOR_EACH_I( PRINT_SEQ_detail, ~, SEQ) std::cout << "\n"; }while(0)
#define PRINT(X) PRINT_SEQ((X))
#define PRINT_M(M) do{ std::cout << #M " = \n"; std::cout << M; }while(0)

namespace numerics{

static const int infinity = -1;

using float_t = long double;

template<class T>
struct basic_matrix{
        using value_type = T;
        enum identity_e{ identity };
        basic_matrix(basic_matrix const&)=default;
        basic_matrix& operator=(basic_matrix const&)=default;
        explicit basic_matrix(size_t height = 2, size_t width = 1):
                height_{height},
                width_{width}, 
                data_(width_ * height_ )
        {}
        explicit basic_matrix(size_t height, size_t width, identity_e):
                height_{height},
                width_{width}, 
                data_(width_ * height_ )
        {
                for(size_t i=0;i!=width;++i)
                        (*this)(i,i) = 1;
        }
        friend std::ostream& operator<<(std::ostream& ostr, basic_matrix<T> const& self){
                std::vector<std::vector<std::string> > buf;
                for(size_t i{0};i!=self.height();++i){
                        buf.emplace_back();
                        for(size_t j{0};j!=self.width();++j){
                                buf.back().emplace_back(boost::lexical_cast<std::string>(self(i,j)));
                        }
                }
                std::vector<size_t> widths(self.width(),0);
                for( auto const& line : buf){
                        for(size_t i=0;i!=self.width();++i){
                                widths[i] = std::max(widths[i], line[i].size());
                        }
                }

                for( auto const& line : buf){
                        ostr << "|";
                        for(size_t i=0;i!=self.width();++i){
                                size_t pad{ widths[i] - line[i].size() };
                                size_t left{ pad / 2};
                                size_t right{ pad - left };
                                ostr 
                                        << ( i == 0 ? "" : " " ) 
                                        << std::string(left,' ') << line[i] << std::string(right,' ');

                        }
                        ostr << "|\n";
                }
                return ostr;
        }
        T const& operator()(size_t i, size_t j)const{
                return data_[i * width_ + j];
        }
        T& operator()(size_t i, size_t j){
                //PRINT_SEQ((i)(j));
                return data_[i * width_ + j];
        }
        T const& operator()(size_t i)const{
                return data_[i];
        }
        T& operator()(size_t i){
                return data_[i];
        }
        size_t width()const{ return width_; }
        size_t height()const{ return height_; }
        T det()const{
                using std::get;
                T ret{};        
                std::vector< std::tuple<size_t, int, T> > aux;
                enum{ E_N, E_Mask, E_Sigma};
                aux.emplace_back(0,0,1);
                for(;aux.size();){
                        auto top{ aux.back() };
                        for( auto const& t : aux){
                                //PRINT_SEQ((get<E_N>(top))(get<E_Mask>(top))(get<E_Sigma>(top)));
                        }
                        aux.pop_back();
                        if( get<E_N>(top) == width_ ){
                                //PRINT_SEQ((get<E_N>(top))(get<E_Mask>(top))(get<E_Sigma>(top)));
                                ret += get<E_Sigma>(top);
                                continue;
                        }
                        int mask = 1;
                        T sign{ +1 };
                        for(size_t i=0;i<width_;++i, mask<<=1){
                                if( !(get<E_Mask>(top) & mask)){
                                        aux.emplace_back( get<E_N>(top)+1, 
                                                          get<E_Mask>(top) | mask, 
                                                          get<E_Sigma>(top) * (*this)(i, get<E_N>(top)) * sign);
                                        sign *= -1;
                                }
                        }
                }
                return ret;
        }
        T cof(size_t i, size_t j)const{
                basic_matrix<T> aux(height()-1, width()-1);
                for(size_t x=0, mi=0;x!=height();++x){
                        if( i == x )
                                continue;
                        for(size_t y=0,mj=0;y!=width();++y){
                                if( j == y )
                                        continue;
                                aux(mi,mj) = (*this)(x,y);
                                ++mj;
                        }
                        ++mi;
                }
                return aux.det();
        }
        basic_matrix<T> adj()const{
                return cof_matrix().transpose();
        }
        basic_matrix<T> cof_matrix()const{
                basic_matrix<T> aux(height(), width());
                for(size_t i{0};i!=height();++i){
                        for(size_t j{0};j!=width();++j){
                                aux(i,j) = this->cof(i,j)  * std::pow(-1, i + j + 2 );
                        }
                }
                return aux;
        }
        basic_matrix<T> inverse()const{
                auto tmp{ adj() };
                tmp /= det();
                return tmp;
        }
        basic_matrix<T> transpose()const{
                basic_matrix<T> aux(width(), height());
                for(size_t i{0};i!=height();++i){
                        for(size_t j{0};j!=width();++j){
                                aux(j,i) = (*this)(i,j);
                        }
                }
                return aux;
        }

        basic_matrix<T>& operator/=(T const& val){
                for( auto& _ : data_)
                        _ /= val;
                return *this;
        }
        friend
        basic_matrix<T> operator*(basic_matrix<T> const& left, basic_matrix<T> const& right){
                assert( left.width() == right.height() && "precondtion failed");
                basic_matrix<T> ret( left.height(), right.width() );
                for(size_t i{0};i!=ret.height();++i){
                        for(size_t j{0};j!=ret.width();++j){
                                for(size_t k=0;k!=left.width();++k){
                                        ret(i,j) += left(i,k) * right(k,j);
                                }
                        }
                }
                return ret;
        }
        basic_matrix<T>& operator-=(basic_matrix<T> const& that){
                assert( width() == that.width() && "precondtion failed");
                assert( height() == that.height() && "precondtion failed");
                for(size_t i{0};i!=height();++i){
                        for(size_t j{0};j!=width();++j){
                                (*this)(i,j) -= that(i,j);
                        }
                }
                return *this;
        }
        friend
        basic_matrix<T> operator-(basic_matrix<T> left, basic_matrix<T> const& right){
                left -= right;
                return std::move(left);
        }
        bool is_vector()const{
                return width() == 1;
        }
        template<int P>
        T norm()const{
                using std::pow;
                if( is_vector() ){
                        if( P == infinity ){
                                return std::fabs(*std::max_element( data_.begin(), data_.end(), [](auto const& l, auto const& r){
                                        return std::fabs(l) < std::fabs(r);
                                }));
                        } else {
                                T ret = T();
                                for(size_t i{0};i!=height();++i){
                                        ret += pow( std::fabs((*this)(i,0)), static_cast<T>(P) );
                                }
                                ret = pow( ret, 1 / static_cast<T>(P) );
                                return ret;
                        }
                } else{
                        // not implemented
                        return -1;
                }
        }
private:
        size_t height_;
        size_t width_;
        std::vector<T> data_;
};


template<class Matrix_Type>
void LU_factorize(Matrix_Type const& a, Matrix_Type& L, Matrix_Type& U){
        L = Matrix_Type(a.height(), a.width(), Matrix_Type::identity );
        U = Matrix_Type(a.height(), a.width());


        enum Op{ Op_L, Op_U};
        std::vector<std::tuple<Op, int, int> > sch;
        for(size_t k=0;k!=a.width();++k){

                for(size_t j=k;j<a.width();++j){
                        sch.emplace_back( Op_U, k, j);
                }

                for(size_t i=k+1;i<a.height();++i){
                        sch.emplace_back( Op_L, i, k);
                }
        }
        using std::get;
        for( auto const& t : sch ){
                auto i = get<1>(t);
                auto j = get<2>(t);
                if( get<0>(t) == Op_U ){
                        //PRINT_SEQ(('U')(i)(j));
                        U(i,j) = a(i,j);
                        for(size_t k = 0; k  < i; ++k)
                                U(i,j) -= L(i,k) * U(k,j);
                } else{
                        //PRINT_SEQ(('L')(i)(j));
                        L(i,j) += a(i,j);
                        for(size_t k = 0; k  < j; ++k)
                                L(i,j) -= L(i,k) * U(k,j);
                        L(i,j) /= U(j,j);
                }
        }
}

/*
        A x = B , x \in R^n
*/
template<class Matrix_Type>
Matrix_Type linear_solve(Matrix_Type const& A, Matrix_Type const& B){
        using value_type = typename Matrix_Type::value_type;
        /*
                A x = B
                L U x = B

                L y = B
                U x = y;
         */
        Matrix_Type L,U,Y(A.height()),X(A.height());
        LU_factorize(A,L,U);


        for(size_t i{0};i!=Y.height();++i){
                value_type val{ B(i) };
                for(size_t j{0};j<i;++j)
                        val -= Y(j) * L(i,j);
                Y(i) = val / L(i,i);
        }

        for(size_t i{X.height()};i!=0;){
                --i;
                value_type val{ Y(i) };
                for(size_t j{X.height()};j!=0 && j > i;){
                        --j;
                        val -= X(j) * U(i, j);
                }
                X(i) = val / U(i,i);
        }

        return X;
}

void test0(){
        basic_matrix<float_t> A(3,3);
        basic_matrix<float_t> I(3,3, basic_matrix<float_t>::identity );
        basic_matrix<float_t> I2(I);
        I2(0,0) = 0;
        I2(1,1) = 0;
        I2(0,1) = 1;
        I2(1,0) = 1;
        basic_matrix<float_t> B(3,1);
        A(0,0) = 1;
        A(0,1) = 1;
        A(0,2) = 1;
        A(1,0) = 2;
        A(1,1) = 4;
        A(1,2) = 2;
        A(2,0) = -1;
        A(2,1) = 5;
        A(2,2) = -4;


        B(0) = 6;
        B(1) = 16;
        B(2) = -3;

        PRINT_M(A);
        PRINT_M(B);

        auto x = linear_solve(A,B);
        PRINT( A* x );
}

template<class X_t, class F_t, class J_t>
auto newton_method_solve(X_t& x0, F_t const& F, J_t const& J){
        /*
                x_{n+1} = x_{n} - J(x_n)^-1 * F(x)
         */
        X_t x{x0};
        for(;;){
                X_t next{ x - J(x).inverse() * F(x) };
                auto d( x - next);
                auto norm{ d.template norm<infinity>() };
                x = next;
                if( norm < 1e-15){
                        break;
                }
        }
        return std::move(x);

}

template<class X_t, class G_t>
auto solve_simultaneous_iteration(X_t& x0, G_t const& G){
        /*
                x_{n+1} = F(x_n)
         */
        X_t x{x0};
        for(;;){
                X_t next{G(x)};
                auto d( x - next);
                x = next;
                auto norm{ d.template norm<infinity>() };
                if( norm < 1e-15){
                        break;
                }
        }
        return std::move(x);
}

void test1(){

        using mat_t = basic_matrix<float_t>;
        using fun_t = std::function<float_t(mat_t const&)>;


        basic_matrix<fun_t> G(2);
        G(0) = [](mat_t const& X)->float_t{ return std::pow( 1-std::pow(X(1),2.0), .5); };
        G(1) = [](mat_t const& X)->float_t{ return 1.0 / std::sqrt(21.0) * std::pow(9.0 - 5.0 * std::pow(X(0),2.0), .5); };
        
        basic_matrix<float_t> X0(2);
        X0(0) = .5;
        X0(1) = .3;

        #if 0
        basic_matrix<fun_t> J(2, 2);
        J(0,0) = [](mat_t const& X)->float_t{ return 0; };
        J(1,0) = [](mat_t const& X)->float_t{ return -5/std::sqrt(12) * X(0) * std::pow(9 - 5 * std::pow(X(1),2.0), -0.5); };
        J(0,1) = [](mat_t const& X)->float_t{ return -X(1) * std::pow(1 - std::pow(X(2),2), -.5); };
        J(1,1) = [](mat_t const& X)->float_t{ return 0; };
        #endif


        auto X{X0};

        auto make_aux = [&](auto const& M){
                return [&](auto const& X)->mat_t{
                        mat_t ret{ M.height(), M.width() };
                        for(size_t i{0};i!=M.height();++i){
                                for(size_t j{0};j!=M.width();++j){
                                        ret(i,j) = M(i,j)(X);
                                }
                        }
                        return std::move(ret);
                };

        };
        auto G_aux{ make_aux(G) };
        auto sret = solve_simultaneous_iteration( X, G_aux );
        PRINT_M(sret);

}


void test2(){

        using mat_t = basic_matrix<float_t>;
        using fun_t = std::function<float_t(mat_t const&)>;

        basic_matrix<float_t> X0(2);
        X0(0) = .5;
        X0(1) = .3;

        /*
                
                f_1 =   x_1^2 +    x_2 ^ 2 - 1
                f_2 = 5 x_1^2 + 21 x_2 ^ 2 - 9
         */
        basic_matrix<fun_t> F(2);
        F(0) = [](mat_t const& X)->float_t{ return  std::pow(X(0), 2.0) +    std::pow(X(1), 2.0) - 1; };
        F(1) = [](mat_t const& X)->float_t{ return 5*std::pow(X(0), 2.0) + 21*std::pow(X(1), 2.0) - 9; };
        

        basic_matrix<fun_t> J(2, 2);
        J(0,0) = [](mat_t const& X)->float_t{ return  2 * X(0); };
        J(1,0) = [](mat_t const& X)->float_t{ return 10 * X(0); };
        J(0,1) = [](mat_t const& X)->float_t{ return 2 * X(1); };
        J(1,1) = [](mat_t const& X)->float_t{ return 42 * X(1); };


        auto X{X0};

        auto make_aux = [&](auto const& M){
                return [&](auto const& X)->mat_t{
                        mat_t ret{ M.height(), M.width() };
                        for(size_t i{0};i!=M.height();++i){
                                for(size_t j{0};j!=M.width();++j){
                                        ret(i,j) = M(i,j)(X);
                                }
                        }
                        return std::move(ret);
                };

        };
        auto J_aux{ make_aux(J) };
        auto F_aux{ make_aux(F) };
        auto nret = newton_method_solve( X, F_aux, J_aux );
        PRINT_M(nret);

}

} // numerics




int main(){
        using namespace numerics;
        test1();
        test2();
}
