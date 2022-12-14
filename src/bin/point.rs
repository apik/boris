use clap::Parser;
use serde_derive::Deserialize;
use std::fs;
use std::process::exit;

use std::fmt;

/// Single point mode
#[derive(Parser, Debug)]
#[clap(name = "BORIS")]
#[clap(author = "Andrey Pikelner <pikelner@theor.jinr.ru>")]
// #[clap(version = "1.0")]
#[clap(about = "Resummation with conformal mapping", long_about = None)]
struct Args {
    /// Resummation parameter "a", conformal transformation
    #[clap(short, value_parser)]
    a: Option<f64>,
    /// Resummation parameter "b", shift in Borel transform
    #[clap(short, value_parser)]
    b: Option<f64>,
    /// Resummation parameter "lam", prefactor exponent
    #[clap(short, value_parser)]
    lam: Option<f64>,
    /// Resummation parameter "q", expansion variable transform
    #[clap(short, value_parser)]
    q: Option<f64>,
    /// Configuration file name
    #[clap(value_parser)]
    input: String,

}

#[derive(Deserialize)]
struct ConfigInput {
    info: Option<String>,
    f:    Option<Vec<f32>>,
    ep:   Option<f32>,
    a:    Option<f32>,
    b:    Option<f32>,
    lam:  Option<f32>,
    q:    Option<f32>,
    n:    Option<i32>,
}

// Save info about input parameters
enum InputTypeFrom {
    Config,
    Args,
    Defaults,
}

impl fmt::Display for InputTypeFrom {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            InputTypeFrom::Config   => write!(f, "CNFG"),
            InputTypeFrom::Args     => write!(f, "ARGS"),
            InputTypeFrom::Defaults => write!(f, "DFLT"),
        }
    }
}


struct InputParameters {
    info: (String, InputTypeFrom),
    f:    (Vec<f64>, InputTypeFrom),
    ep:   (f64, InputTypeFrom),
    a:    (f64, InputTypeFrom),
    b:    (f64, InputTypeFrom),
    lam:  (f64, InputTypeFrom),
    q:    (f64, InputTypeFrom),
}

fn main() {

    // Parse command line arguments
    let args = Args::parse();
    let contents = match fs::read_to_string(args.input.clone()) {
        Ok(c) => c,
        Err(_) => {
            eprintln!("Could not read file `{}`", args.input);
            exit(1);
        }
    };

    let cfg_data: ConfigInput = match toml::from_str(&contents) {
        Ok(cfg) => cfg,
        Err(_) => {
            eprintln!("Unable to load data from `{}`",args.input);
            exit(1);
        }
    };

    // Initialize default values
    let mut ip = InputParameters{
        info: ("".to_string(), InputTypeFrom::Defaults),
        f: (Vec::new(), InputTypeFrom::Config),
        ep: (0f64, InputTypeFrom::Config),
        a: (0f64, InputTypeFrom::Defaults),
        b: (0f64, InputTypeFrom::Defaults),
        lam: (0f64, InputTypeFrom::Defaults),
        q: (0f64, InputTypeFrom::Defaults),
    };

    // 
    // Input from config file only
    // 
    ip.info = (match cfg_data.info {
        Some(p) => p,
        None => "".to_string(),
    }, InputTypeFrom::Config);
    // input series f = [...]
    ip.f = (cfg_data.f.expect("Provide input series in the form f=[0.1, 0.2, 0.3,...]").iter().map(|&val| val as f64).collect(), InputTypeFrom::Config);
    // input ep
    ip.ep = (cfg_data.ep.expect("Provide ep=...") as f64, InputTypeFrom::Config);

    //
    // Input from file and command line
    //

    ip.a = match args.a {
        Some(a) => (a, InputTypeFrom::Args),
        None    => {
            match cfg_data.a {
                Some(a) => (a as f64, InputTypeFrom::Config),
                None    => (0f64, InputTypeFrom::Defaults),
            }
        },
    };

    ip.b = match args.b {
        Some(b) => (b, InputTypeFrom::Args),
        None    => {
            match cfg_data.b {
                Some(b) => (b as f64, InputTypeFrom::Config),
                None    => (0f64, InputTypeFrom::Defaults),
            }
        },
    };

    ip.lam = match args.lam {
        Some(lam) => (lam, InputTypeFrom::Args),
        None    => {
            match cfg_data.lam {
                Some(lam) => (lam as f64, InputTypeFrom::Config),
                None    => (0f64, InputTypeFrom::Defaults),
            }
        },
    };

    ip.q = match args.q {
        Some(q) => (q, InputTypeFrom::Args),
        None    => {
            match cfg_data.q {
                Some(q) => (q as f64, InputTypeFrom::Config),
                None    => (0f64, InputTypeFrom::Defaults),
            }
        },
    };    

    // Print summary
    println!("[INFO] {}", ip.info.0);
    println!("       ep   = {}", &ip.ep.0);
    println!("[{}] a   = {}", ip.a.1, &ip.a.0);
    println!("[{}] b   = {}", ip.b.1, &ip.b.0);
    println!("[{}] lam = {}", ip.lam.1, &ip.lam.0);
    println!("[{}] q   = {}", ip.q.1, &ip.q.0);
    println!("Input series:");
    for (i, t) in ip.f.0.iter().enumerate() {
        println!("  x^{ord}: {num:.prec$}", prec = 5, num = t, ord = i);
    }

    let integ_ord = cfg_data.n.unwrap_or(20);



    let rsm:f64  = boris::resummator::Resummator::new
        (integ_ord as usize,
         ip.a.0,       // a
         ip.b.0,       // b
         ip.ep.0,      // ep
         ip.f.0.as_mut_slice()
        ).point(ip.lam.0, ip.q.0);

    // Output result
    println!("res = {}", rsm);

}


#[cfg(test)]
#[macro_use]
extern crate assert_float_eq;
mod tests {
    // ETA
    #[test]
    fn test_point_eta_0() {

        let rsm = boris::resummator::Resummator::new
            (20,
             3f64/8f64, // a
             11.0,      // b
             1f64,      // ep
             &[0.0, 0.0, 0.015625, 0.0166015625, -0.008367494045052055, 0.026504617486672634, -0.090730434012201]).point(2.54, 0.10);
        assert_float_relative_eq!(rsm, 0.03100, 0.001); 

    }

    #[test]
    fn test_point_eta_1() {

        let rsm = boris::resummator::Resummator::new
            (20,
             1f64/3f64, // a
             11.0,      // b
             1f64,      // ep
             &[0.0, 0.0, 0.018518518518518517, 0.018689986282578876, -0.00832877034109173, 0.025656450831775084, -0.08127264250445679]).point(2.56, 0.10);
        assert_float_relative_eq!(rsm, 0.03615, 0.001);

    }

    #[test]
    fn test_point_eta_2() {

        let rsm = boris::resummator::Resummator::new
            (20,
             3f64/10f64, // a
             11.5,       // b
             1f64,       // ep
             &[0.0, 0.0, 0.02, 0.019, -0.007893594032531367, 0.023209119705921755, -0.06862674086716794]).point(2.52, 0.10);
        assert_float_relative_eq!(rsm, 0.03798, 0.001); 

    }

    #[test]
    fn test_point_eta_3() {

        let rsm = boris::resummator::Resummator::new
            (20,
             3f64/11f64, // a
             12.0,       // b
             1f64,       // ep
             &[0.0, 0.0, 0.02066115702479339, 0.01839867495389659, -0.007449474404300666, 0.02038282090515603, -0.057024194932236365]).point(2.56, 0.09);
        assert_float_relative_eq!(rsm, 0.03783, 0.001);

    }

    #[test]
    fn test_point_eta_4() {

        let rsm = boris::resummator::Resummator::new
            (20,
             1f64/4f64, // a
             14.0,      // b
             1f64,      // ep
             &[0.0, 0.0, 0.020833333333333332, 0.017361111111111112, -0.007085182272194144, 0.017631497114137844, -0.047362839351141696]).point(2.52, 0.08);
        assert_float_relative_eq!(rsm, 0.03663, 0.001);

    }

}
