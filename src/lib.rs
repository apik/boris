
// Integrator
pub mod integrator;

// Resummator
pub mod resummator;



#[cfg(test)]
#[macro_use]
extern crate assert_float_eq;

mod tests {

    #[test]
    fn test_gen_laguerre() {
        // |(a - b) / a| > epsilon
        assert_float_relative_eq!(crate::integrator::gen_laguerre(3, 1.3 , 10.0), -17.1772, 0.001);
        assert_float_relative_eq!(crate::integrator::gen_laguerre(3, 0.0 , 7.0), -3.66667, 0.001);
    }


    #[test]
    fn test_imtqlx() {
        let mut d: Vec<f64> = vec![1.0, 2.0, 3.0];
        let mut e: Vec<f64> = vec![0.3, 0.1, 0.0];
        let mut z: Vec<f64> = vec![0.0, 1.0, 0.0];
        assert!(crate::integrator::imtqlx(3, d.as_mut_slice(), e.as_mut_slice(), z.as_mut_slice()).is_ok());
        // d
        assert_float_relative_eq!(d[0], 0.916561, 0.001);
        assert_float_relative_eq!(d[1], 2.07308, 0.001);
        assert_float_relative_eq!(d[2], 3.01036, 0.001);
        // z
        assert_float_relative_eq!(z[0], -0.267935, 0.001);
        assert_float_relative_eq!(z[1], 0.957915, 0.001);
        assert_float_relative_eq!(z[2], 0.103001, 0.001);
    }


    #[test]
    fn test_gauss_laguerre_weights() {
        let mut weights: Vec<f64>   = vec![0f64;5];
        let mut nodes: Vec<f64>     = vec![0f64;5];
        assert!(crate::integrator::gauss_laguerre_weights(5,1.3, nodes.as_mut_slice(), weights.as_mut_slice()).is_ok());
        assert_float_relative_eq!(weights[4], 0.000100544, 0.001);
        assert_float_relative_eq!(nodes[4], 14.7354, 0.001);
    }


    #[test]
    fn test_gauss_laguerre_struct() {
        // n=20, alpha = 4.5
        let gl = crate::integrator::GaussLaguerre::new(20, 4.5);
        // check only last element
        assert_float_relative_eq!(gl.w(19), 1.29619e-23, 0.001);
        assert_float_relative_eq!(gl.x(19), 74.7087, 0.001);
    }


    #[test]
    fn test_poch_div_fac() {
        assert_float_relative_eq!(crate::resummator::poch_div_fac(0.1, 5), 0.0244668, 0.001);
        assert_float_relative_eq!(crate::resummator::poch_div_fac(4.3, 0), 1.0, 0.001);
        assert_float_relative_eq!(crate::resummator::poch_div_fac(4.3, 1), 4.3, 0.001);
        assert_float_relative_eq!(crate::resummator::poch_div_fac(4.3, 2), 11.395, 0.001);
        assert_float_relative_eq!(crate::resummator::poch_div_fac(4.3, 7), 165.339, 0.001);
    }

}
