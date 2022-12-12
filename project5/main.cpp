// FYS4150
//
// Project 5
//
// Solving Schroedinger equation

# include <iostream>
# include <armadillo>
# include <complex>
# include <string>

using namespace std;
using namespace arma;

complex<double> i(0.,1.); // imaginary number

void normalise_state(cx_mat& U)
{
    // normalises state matrix U using the sum of |u_ij|^2
    double s = real(accu(U % conj(U)));
    cout << "s = " << s << endl;
    U = U/sqrt(s);
}

cx_mat initialise_wave_packet( int n_inner, double h,
    double x_center, double x_spread, double p_x,
    double y_center, double y_spread, double p_y)
{
    // initialises state using a Gaussian wave packet
    cx_mat U(n_inner, n_inner, fill::zeros);

    // declare variables used in Gaussian function
    double x, y;
    double x_term, y_term;
    complex<double> px_term, py_term;

    // loop over x
    for (int i_x = 0; i_x < n_inner; i_x ++)
    {
        x = i_x*h;
        x_term = -0.5*(x-x_center)*(x-x_center)/(x_spread*x_spread);
        px_term = i*p_x*(x-x_center);

        // loop over y
        for (int i_y = 0; i_y < n_inner; i_y ++)
        {
            y = i_y*h;
            y_term = -0.5*(y-y_center)*(y-y_center)/(y_spread*y_spread);
            py_term = i*p_y*(y-y_center);
            U(i_x,i_y) = exp(x_term + y_term + px_term + py_term);
        }
    }
    
    normalise_state(U); // normalise state
    
    return U;
}

cx_mat create_A_matrix(int n_inner, cx_vec a, complex<double> r)
{
    cx_mat A(n_inner*n_inner, n_inner*n_inner);
    //A.diag() = a;
    //A.diag(n_inner).fill(-r);
    int sub_idx, diag_idx;
    for (sub_idx=0; sub_idx<n_inner; sub_idx++)
    {
        // looping over submatrices

        for (diag_idx=0; diag_idx < n_inner - 1; diag_idx++)
        {
            // looping over diagonal in submatrix
            // set diagonal A(i,i)
            A(sub_idx*n_inner + diag_idx, sub_idx*n_inner + diag_idx) = a(sub_idx*n_inner + diag_idx);
            // set superdiagonal of submatrix A(i,i+1)
            A(sub_idx*n_inner + diag_idx, sub_idx*n_inner + diag_idx + 1) = -r;
            if (sub_idx < n_inner - 1)
            {
                // set superduperdiagonal A(i,i+M-2) (not valid for final submatrix, hence an if test)
                A(sub_idx*n_inner + diag_idx, (sub_idx+1)*n_inner) = -r;
            }
        }
        // set final diagonal element in submatrix
        A(sub_idx*n_inner +  diag_idx, sub_idx*n_inner + diag_idx) = a(sub_idx*n_inner + diag_idx);
        
    }
    A = symmatu(A, false); // mirroring upper and lower triangular part
    return A;
}

int main(int argc, char* argv[])
{
    string potential = argv[1];
    double h = stod(argv[2]);
    double dt = stod(argv[3]);
    double T = stod(argv[4]);
    double x_center = stod(argv[5]);
    double x_spread = stod(argv[6]);
    double p_x = stod(argv[7]);
    double y_center = stod(argv[8]);
    double y_spread = stod(argv[9]);
    double p_y = stod(argv[10]);

    // loading potential matrix
    mat V;
    V.load(potential, raw_ascii);

    int n_inner = V.n_cols; // number of inner points = M-2
    h = 1./(n_inner+1);

    cout << "h = " << h << endl;

    // creating A and B matrices
    complex<double> r = i*dt/(2.*h*h);
    cx_vec a = ones(n_inner*n_inner)*(1. + 4.*r) + 0.5*i*dt*V.as_col();
    cx_mat A = create_A_matrix(n_inner, a, r);
    cx_vec b = ones(n_inner*n_inner)*(1. - 4.*r) - 0.5*i*dt*V.as_col();
    cx_mat B = create_A_matrix(n_inner, b, -r);

    cx_mat U;
    U = initialise_wave_packet(n_inner, h, x_center, x_spread, p_x, y_center, y_spread, p_y);

    cx_vec u, u_;
    u = U.as_col();
    cx_vec Bu;

    int n_timepoints = T/dt + 1;
    
    cube U_squared = cube(n_inner, n_inner, n_timepoints);
    cube U_real = cube(n_inner, n_inner, n_timepoints);
    cube U_imag = cube(n_inner, n_inner, n_timepoints);
    cout << "u size " << size(u) << endl;
    cout << n_inner*n_inner << endl;
    cout << "DEBUG 1" << endl;
    //U = reshape()
    U_squared.slice(0) = real( U % conj(U) );
    U_real.slice(0) = real(U);
    U_imag.slice(0) = imag(U);

    for (int t=1; t<n_timepoints; t++)
    {
        Bu = B*u;
        u_ = solve(A, Bu);
        u = u_;
        U = reshape(u, n_inner, n_inner);
        cout << 100.*t/n_timepoints << "%" << endl;

        U_squared.slice(t) = real( U % conj(U) );
        U_real.slice(t) = real(U);
        U_imag.slice(t) = imag(U);
    }
    

    U_squared.save("p.bin");
    U_real.save("real.bin");
    U_imag.save("imag.bin");
    //cout << A << endl;
    // mat A_ = ones(n_inner*n_inner, 4);
    // A_(3,2) = 10;
    // cout << a << endl;
    // cout << endl;
    // vec A_vec = A_.as_col();
    // cout<< A_vec << endl;
    // cout << A_vec(3 + 2*n_inner*n_inner) << endl;

    return 0;
}

/*

read input file

set up potential matrix: single, double, triple slit

set up A, B

solve Ax = Bx_

*/

// int i, j, k;
// int n_inner;
// k = i + j*n_inner;
// u11 u12 u13
// u21 u22 u23
// u31 u32 u33

// u11 u21 u31 --- u23 u33


// cx_vec a = ones(n_inner*n_inner)*(1. + 4.*r) + 0.5*i*dt*V.as_col();
// sp_cx_mat A = create_A_matrix(n_inner, a, r);
// cx_vec b = ones(n_inner*n_inner)*(1. - 4.*r) - 0.5*i*dt*V.as_col();
// sp_cx_mat B = create_A_matrix(n_inner, b, -r);

// cx_mat U = initialise_wave_packet(n_inner, h, x_center, x_spread, p_x, y_center, y_spread, p_y);
// cx_vec u = U.as_col();

// int n_timepoints = T/dt + 1;

// cube U_squared = cube(n_inner, n_inner, n_timepoints);
// cube U_real = cube(n_inner, n_inner, n_timepoints);
// cube U_imag = cube(n_inner, n_inner, n_timepoints);

// U_squared.slice(0) = real(U % conj(U));
// U_real.slice(0) = real(U);
// U_imag.slice(0) = imag(U);

// cx_vec Bu;

// for (int t = 1; t<n_timepoints; t++)
// {
//     Bu = B*u;
//     spsolve(u, A, Bu);
//     cout << "finished solve" << endl;
//     U = reshape(u, n_inner, n_inner);
//     U_squared.slice(t) = real(U % conj(U));
//     U_real.slice(t) = real(U);
//     U_imag.slice(t) = imag(U);
// }