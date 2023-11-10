// BGW using Shamir secret on linear circuit
use rand::prelude::*;

fn gen_ran_in_fp(prime:i64)-> i64{
    let mut rng = rand::thread_rng(); //why thread_rng
    let y:i64 = rng.gen();
    (prime+y%prime)%(prime)
}

fn gen_poly(m:i64, prime:i64 , degree:u8) -> Vec<i64>{
    let mut vec:Vec<i64> = Vec::with_capacity((degree+1) as usize);
    vec.push(m);
    for _ in 0..(degree) {
        //let tmp = gen_ran_in_fp(prime);
        vec.push(gen_ran_in_fp(prime));
    } 
    vec
}

#[allow(dead_code)]
fn eval_px_at_xi(prime:i64, pn:&Vec<i64>, x:i64)-> i64{
    let mut y:i64 = (*pn)[0]; //init
    let mut x_pow:i64 = x;

    for i in 1..(*pn).len() {
        let tmp = (x_pow*(*pn)[i])%prime;

        y = (y + tmp)%prime;
        x_pow = (x_pow*x)%prime;
        //println!("{} {}",xi,x_pow);
    }
    y
}   

#[allow(dead_code)]
fn gen_share(prime:i64, pn:&Vec<i64>, xi:&Vec<i64>, no_of_shares:i64) -> Vec<i64>{ //one share at a time
    let mut yi:Vec<i64> = Vec::with_capacity(no_of_shares as usize);
    
    // evaluating at 1,2,3,4 when Field over 5
    for i in 0..(no_of_shares as usize) {
        yi.push(eval_px_at_xi(prime, &pn, xi[i]));
    }
    yi
}

fn gcd(x:i64,y:i64) -> i64{
    if y==0 {
        x
    }else{
        gcd(y,x%y)
    }
}
fn inv_modp(a: i64,prime:i64) -> i64{ //change a to take -ve valeus, else sanitize inputs before
    if gcd(a,prime) != 1 {
        println!("Not coprime. Can't find inverse.");
        0
    }else { //without else it doesn't work
        let a = (prime+a)%prime; //defined antoher a
        let mut r:[i64;3] = [a,prime,0];
        let mut x:[i64;3] = [1,0,0];
        let mut y:[i64;3] = [0,1,0];
        //let mut q:i64=0;
        
        while r[1] != 1i64 {
            r[2] = (r[0] % r[1]) % prime;
            let q = (r[0] / r[1]) % prime;
            x[2] = (x[0] + (prime- q)*x[1]) % prime; // x[0] - q*x[1] in i64
            y[2] = (y[0] + (prime- q)*y[1]) % prime;

            r[0] = r[1]; r[1] = r[2];
            x[0] = x[1]; x[1] = x[2];
            y[0] = y[1]; y[1] = y[2];

        }
        x[1]
    }
}

#[allow(dead_code)]
fn reconstruct(xi:&Vec<i64>, yi:&Vec<i64>,prime:i64) -> i64{
    /* Steps
        1. Lagrange interpolation
        2. Return m, the constant coefficient
     */
    let n = xi.len();
    let mut res:i64 = 0;

    for i in 0..n {
        let mut num:i64=1;
        let mut den:i64=1;
        for j in 0..n {
            if i == j {continue;}
            else {
                let tmpj:i64 = prime-xi[j];  
                num = (num * (tmpj)) % prime;
                den = (den * ((xi[i] + tmpj) % prime)) % prime; 
            }
        }
        res = (res + (yi[i] * (num * inv_modp(den, prime))%prime)%prime)%prime;
    }

    res
}

/* N parties
 * Party 1 has pvtInput A, Party 2 has B, C, D accordingly.
 * CIrcuit -> K = A,B // non linear
 * Steps:
 * 1. Secret share A, B to get Ai's and Bi's
 * 2. Multiply coefficent wise to get Ki's = Ai * Bi
 * 3. Secret share each Ki so that AB is a linear fn of Ki's; AB = c1.K1 + c2.K2 + .. cn.Kn where ci's are const
 * 4. Now apply (n,t) shamir secret for each Ki
 * 5.  
 */

fn main(){
    const prime:i64=5;
    const no_of_parties:i64 = 4;
    const no_of_shares:i64=4;
    const degree:u8 = 2;

    // f(A,B,C,D) = A.B + C.D
    let mut pvt_inputs: Vec<i64> = Vec::with_capacity(no_of_parties as usize); // A, B, C, D
    pvt_inputs.push(4); pvt_inputs.push(2); pvt_inputs.push(0);pvt_inputs.push(3); 

    let mut xi:Vec<i64> = Vec::with_capacity(no_of_shares as usize);
    for i in 0..(no_of_shares as usize) {
        /* need to check for unique x's*/
        let mut tmp:i64 = gen_ran_in_fp(prime); //gen xi
        while tmp == 0 || xi.contains(&tmp) { // linear op on all Fp
            tmp = gen_ran_in_fp(prime);
        } 
        xi.push(tmp);
    }
    //print!("{} {} {} {}",xi[0],xi[1],xi[2],xi[3]);
    let mut shares:[Vec<i64>;no_of_parties as usize] = Default::default();
    

    //STEP 1 -> GEN SHARES for its own pvt_input
    for i in 0..(no_of_parties as usize) {
        let pn:Vec<i64> = gen_poly(pvt_inputs[i], prime, degree);
    
        for _ in 0..(no_of_shares as usize) {
           shares[i as usize] = gen_share(prime, &pn, &xi, no_of_shares);
        }
        //print!("{} {} {} {}\n",shares[i][0],shares[i][1],shares[i][2],shares[i][3]);
    }

    //STEP 2 -> share the generated shares
    //STEP 3 -> add the values to get partial yi's
    let mut rshares:[Vec<i64>;no_of_parties as usize] = Default::default();

    let mut partial_yi:Vec<i64> = Vec::with_capacity(no_of_parties as usize);
    for i in 0..(no_of_parties as usize){
        let mut tmp_sum:i64=0;
        for j in 0..(no_of_shares as usize){
            rshares[i].push(shares[j][i]);
            tmp_sum = (tmp_sum + shares[j][i])%prime;
        }
        partial_yi.push(tmp_sum);
        print!("{} {}\n",xi[i],partial_yi[i]);
        //print!("{} {} {} {}\n",rshares[i][0],rshares[i][1],rshares[i][2],rshares[i][3]);
    }

    //STEP 4 -> Recon using rshares to get Y
    let mdash = reconstruct(&xi,&partial_yi, prime);
    println!("{}",mdash);
    
}
