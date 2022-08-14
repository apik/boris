use special::Gamma;
// Matrix diagonalisation
pub fn imtqlx ( n: usize, d: &mut [f64], e: &mut [f64], z: &mut [f64]) -> Result<String, String> {

    const MAXITER:usize = 30;
    return match n {

        1 => Ok("n=1".to_string()),

        _ => {
            e[n-1] = 0.0;
            
            let mut m = 0;
            
            for l in 1.. (n+1)  {
                let mut j = 0;
                loop {
                    for mm in l..(n+1){
                        m = mm;
                        if m == n { break; }
                        if e[m-1].abs() <= f64::EPSILON * ( d[m-1].abs() + d[m].abs() )  { break; }
                    }
                    let mut p = d[l-1];
                    if m == l { break; }
                    if j >= MAXITER { return Err("Maximum iteration exceeded!".to_string()) }
                    j +=  1;
                    let mut g = ( d[l] - p ) / ( 2f64 * e[l-1] );
                    let mut r =  ( g * g + 1f64 ).sqrt();
                    g = d[m-1] - p + e[l-1] / ( g + r.abs().copysign(g) );
                    let mut s = 1f64;
                    let mut c = 1f64;
                    p = 0f64;
                    let mml: usize = m - l;
                    
                    for ii in 1..(mml+1) {
                        let i = m - ii;
                        let mut f = s * e[i-1];
                        let b = c * e[i-1];
                        
                        if g.abs() <= f.abs() {
                            c = g / f;
                            r = ( c * c + 1f64 ).sqrt();
                            e[i] = f * r;
                            s = 1f64 / r;
                            c = c * s;
                        } else{
                            s = f / g;
                            r = ( s * s + 1f64 ).sqrt();
                            e[i] = g * r;
                            c = 1f64 / r;
                            s = s * c;
                        }
                        g = d[i] - p;
                        r = ( d[i-1] - g ) * s + 2f64 * c * b;
                        p = s * r;
                        d[i] = g + p;
                        g = c * r - b;
                        f = z[i];
                        z[i] = s * z[i-1] + c * f;
                        z[i-1] = c * z[i-1] - s * f;
                    }
                    d[l-1] = d[l-1] - p;
                    e[l-1] = g;
                    e[m-1] = 0f64;
                }
            }
            /* Sorting. */
            for ii in 2..(m+1) {
                let i = ii - 1;
                let mut k = i;
                let mut p = d[i-1];
                
                for j in ii..(n+1) {
                    if d[j-1] < p {
                        k = j;
                        p = d[j-1];
                    }
                }
                
                if k != i {
                    d[k-1] = d[i-1];
                    d[i-1] = p;
                    p = z[i-1];
                    z[i-1] = z[k-1];
                    z[k-1] = p;
                }
            }
            
            Ok("Success!".to_string())            
        },
    }
}


// LaguerreL[n,a,x]
pub fn gen_laguerre(n: usize, alpha: f64, x: f64) -> f64 {
    match n {
        0 => 1f64,
        1 => 1f64 + alpha - x,
        _ => ((alpha - x - 1f64 + (2*n) as f64)*gen_laguerre(n-1, alpha, x)
              - (alpha - 1f64 + n as f64)*gen_laguerre(n-2, alpha, x) )/(n as f64),
    }
}



// Constructor of Gauss-Laguerre weights for fixed (n,alpha)
// integrating ord(2n-1) polynoimal expression exactly
pub fn gauss_laguerre_weights(n: usize, alpha: f64, nodes: &mut [f64], weights: &mut [f64])  -> Result<String, String> {

    let mut subdiag: Vec<f64>   = vec![0f64;n];

    weights[0] = (alpha+1f64).gamma().sqrt();

    // Fill diagonal and subdiagonal elements of the matrix
    for i in 0..n {
        nodes[i]   = alpha + (2*i + 1) as f64;
        subdiag[i] = ((alpha + (i+1) as f64)*(i+1) as f64).sqrt();
    }

    if imtqlx(n, nodes, subdiag.as_mut_slice(), weights).is_ok() {
        // square weights
        for w in weights.iter_mut() { *w *= *w; }
        return Ok("Success!".to_string())
    }
    return Err("Unable to generate nodes and weights!".to_string())
}

// Caclualate node coordinates and weights
pub struct GaussLaguerre {
    pub n: usize, 
    pub alpha: f64,
    pub nodes: Vec<f64>,
    pub weights: Vec<f64>, 
}

impl GaussLaguerre {
    pub fn new(n: usize, alpha: f64) -> GaussLaguerre {
        let mut weights: Vec<f64>   = vec![0f64; n];
        let mut nodes: Vec<f64>     = vec![0f64; n];
        let _ok = gauss_laguerre_weights(n, alpha, nodes.as_mut_slice(), weights.as_mut_slice());
        GaussLaguerre {
            n,
            alpha,
            nodes,
            weights,
        }
    }
    // Get i-th node
    pub fn x(&self, i: usize) -> f64 { self.nodes[i]   }
    // Get i-th weight
    pub fn w(&self, i: usize) -> f64 { self.weights[i] }
}



