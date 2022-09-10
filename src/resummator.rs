use special::Gamma;
extern crate num;
use num::integer::binomial;
use crate::integrator;

// (a)_n/n!
pub fn poch_div_fac(a: f64, n: usize) -> f64 {

    match n {
        0 => 1_f64,
        1 => a,
        _ => poch_div_fac(a, n-1)*(a + (n-1) as f64)/(n as f64)
    }
}


// Conformal mapping with parameter a
fn cmap(x: f64, a: f64) -> f64 {
    ((1f64+a*x).sqrt() - 1f64)/((1f64+a*x).sqrt() + 1f64)
}




pub struct Resummator<'slice> {
    // Integrateion rule order
    pub n: usize,
    // Fixed, since never changed
    pub a: f64,
    // Fixed, since needs rules recalculation
    pub b: f64,
    // Fixed, since known from input
    pub ep: f64,
    // Series coefficients
    pub f: &'slice [f64],

    pub nodes: Vec<f64>,
    pub weights: Vec<f64>, 
}

impl<'slice> Resummator<'slice> {
    pub fn new(n: usize, a: f64, b: f64, ep: f64, f: &[f64]) -> Resummator {
        let mut weights: Vec<f64>   = vec![0f64; n];
        let mut nodes: Vec<f64>     = vec![0f64; n];
        let _ok = integrator::gauss_laguerre_weights(n, b, nodes.as_mut_slice(), weights.as_mut_slice());

        Resummator {
            n,a,b,ep,f,
            nodes,
            weights,
        }
    }

    pub fn point(&self, lam: f64, q: f64) -> f64 {
        let len: usize = self.f.len();
        // transformed series
        let mut f_t: Vec<f64>   = vec![0f64; len];
        f_t[0] = self.f[0];
        for r in 1..len {
            for k in 0..(r+1) {
                f_t[r] += self.f[k]*(binomial(r-1,r-k) as f64)*((-q).powi((r-k) as i32));
            }
        }


        // New variable

        let ep_t = self.ep/(1f64 - q*self.ep);

        let mut bv: Vec<f64>   = vec![0f64; len];

        // Fill B[r]
        for r in 0..len {
            for k in 0..(r+1) {
                bv[r] += f_t[k]/((self.b + 1f64 + k as f64).gamma())*
                    (4f64.powf(-lam))*((4f64/self.a).powi(k as i32))*
                    (-1f64).powi((r-k) as i32)*
                    poch_div_fac(1f64 + 2f64*lam - (r as f64)  - (k as f64) ,r-k);
            }
        }



        // Function evaluation on the grid
        let mut int: f64 = 0f64;
        // Iterate over nodes
        for k in 0..self.n {
            let t = self.nodes[k];
            let w = cmap(ep_t*t, self.a);
            let pf = (self.a*ep_t*t/w).powf(lam);

            // Function from series in "w"
            let mut f_x: f64 = 0f64;
            for j in 0..len {
                f_x += pf*bv[j]*w.powi(j as i32);
            }
            int += f_x*self.weights[k];
        }
        int
    }
}
