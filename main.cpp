#include <iostream>
#include <vector>
#include <cassert>

#include <boost/lexical_cast.hpp>
#include <boost/preprocessor.hpp>

#define PRINT_SEQ_detail(r, d, i, e) do{ std::cout << ( i ? ", " : "" ) << BOOST_PP_STRINGIZE(e) << " = " << (e); }while(0);
#define PRINT_SEQ(SEQ) do{ BOOST_PP_SEQ_FOR_EACH_I( PRINT_SEQ_detail, ~, SEQ) std::cout << "\n"; }while(0)
#define PRINT(X) PRINT_SEQ((X))
#define PRINT_M(M) do{ std::cout << #M " = \n"; std::cout << M; }while(0)

template<class T>
struct basic_matrix{
        using value_type = T;
        enum identity_e{ identity };
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
        basic_matrix<long double> A(3,3);
        basic_matrix<long double> I(3,3, basic_matrix<long double>::identity );
        basic_matrix<long double> I2(I);
        I2(0,0) = 0;
        I2(1,1) = 0;
        I2(0,1) = 1;
        I2(1,0) = 1;
        basic_matrix<long double> B(3,1);
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


int main(){
        test0();
}
