use rayon::prelude::*;
use clap::Parser;
use serde_derive::Deserialize;
use std::fs;
use std::process::exit;

use std::fs::File;
use std::io::Write;
use indicatif::ProgressBar;
use indicatif::ProgressStyle;
use console::style;

/// Scan mode
#[derive(Parser, Debug)]
#[clap(name = "BORIS - scan mode")]
#[clap(about = None, long_about = None)]
struct Args {
    /// Number of workers for parallel run
    #[clap(short, long, value_parser)]
    workers: Option<usize>,
    /// Resummation tables dump file name
    #[clap(short, long, value_parser)]
    dump: Option<String>,
    /// Output file name
    #[clap(short, long, value_parser)]
    out: Option<String>,
    /// Configuration file name
    #[clap(value_parser)]
    input: String,
}

#[derive(Deserialize)]
struct ConfigInput {
    info:     Option<String>,
    f:        Vec<f32>,
    ep:       f32,
    a:        f32,
    b:        Vec<f32>,
    lam:      Vec<f32>,
    q:        Vec<f32>,
    intord:   Option<i32>,
    dump:     Option<String>,
    out:      Option<String>,
    deltab:   usize,
    deltalam: usize,
    deltaq:   usize,
}

// Triple [l, u]/d
struct Par3 {
    // lower limit
    l: f64,
    // upper limit
    u: f64,
    // divisions
    d: usize,
}

impl Par3 {
    pub fn min(&self) -> f64 {self.l}
    pub fn max(&self) -> f64 {self.u}
    pub fn div(&self) -> usize {self.d}
    pub fn sites(&self) -> usize {self.d + 1}
    pub fn step(&self) -> f64 {(self.u - self.l)/((self.d + 1) as f64)}
}


// Return min from entered index and lower border
fn min_or_cut(pos: usize, delta: usize) -> usize {
    if delta > pos { 0 } else { pos - delta }
}

fn max_or_cut(p: &Par3, pos: usize, delta: usize) -> usize {
    if pos + delta + 1 > p.sites() { p.sites() } else { pos + delta + 1}
}


type Grid3d = Vec<Vec<Vec<f64>>>;
type ErrComponents = (f64, usize, usize, usize, f64, f64, f64, f64, f64);

struct ErrorEstimator {
    // calculated grids
    g0: Grid3d,
    g1: Grid3d,
    g2: Grid3d,
    // parameters
    b: Par3,
    lam: Par3,
    q: Par3,
    // Ranges
    delta_b: usize,
    delta_lam: usize,
    delta_q: usize,
}

impl ErrorEstimator {
    pub fn new(g0: Grid3d, g1: Grid3d, g2: Grid3d,
               b: Par3, lam: Par3, q: Par3) -> ErrorEstimator {
        ErrorEstimator {
            g0, g1, g2, b, lam, q,
            delta_b: 0 , delta_lam: 0, delta_q: 0,
        }
    }

    // 
    // Variation calculation functions
    //


    fn b_var(&self, b_pos: usize, lam_pos: usize, q_pos: usize) -> f64 {
        let f0 = self.g0[b_pos][lam_pos][q_pos];
        let f1 = self.g1[b_pos][lam_pos][q_pos];
        let mut var0: f64 = 0f64;
        let mut var1: f64 = 0f64;
        for ib in min_or_cut(b_pos, self.delta_b)..max_or_cut(&self.b, b_pos, self.delta_b) {
            let df0 = (f0 - self.g0[ib][lam_pos][q_pos]).abs();
            if df0 > var0 {var0 = df0}
            let df1 = (f1 - self.g1[ib][lam_pos][q_pos]).abs();
            if df1 > var1 {var1 = df1}
        }
        // max of two errors
        if var0 > var1 { var0 } else { var1 }
    }

    fn lam_var(&self, b_pos: usize, lam_pos: usize, q_pos: usize) -> f64 {
        let f0 = self.g0[b_pos][lam_pos][q_pos];
        let mut var0: f64 = 0f64;
        for ilam in min_or_cut(lam_pos, self.delta_lam)..max_or_cut(&self.lam, lam_pos, self.delta_lam) {
            let df0 = (f0 - self.g0[b_pos][ilam][q_pos]).abs();
            if df0 > var0 {var0 = df0}
        }
        var0
    }

    fn q_var(&self, b_pos: usize, lam_pos: usize, q_pos: usize) -> f64 {
        let f0 = self.g0[b_pos][lam_pos][q_pos];
        let mut var0: f64 = 0f64;
        for iq in min_or_cut(q_pos, self.delta_q)..max_or_cut(&self.q, q_pos, self.delta_q) {
            let df0 = (f0 - self.g0[b_pos][lam_pos][iq]).abs();
            if df0 > var0 {var0 = df0}
        }
        var0
    }



    pub fn scan(&mut self, delta_b: usize, delta_lam: usize, delta_q: usize) -> (f64, f64, f64) {
        // Store ranges of variation
        self.delta_b = delta_b;
        self.delta_lam = delta_lam;
        self.delta_q = delta_q;

        // err -> (b,lam,q,   f,  d_l, var_b, var_lam, var_q)
        let mut err_db: Vec<ErrComponents> =
            Vec::with_capacity(self.b.sites()*self.lam.sites()*self.q.sites());
        // Calculate error on each site oof the grid
        println!("\n{:^40}", style("SCAN GRID").bold().dim());
        let pb = ProgressBar::new((self.b.sites()*self.lam.sites()*self.q.sites()) as u64);
        pb.set_style(ProgressStyle::default_bar().template("{spinner:.gray} {bar:40}").unwrap());
        for b_pos in 0..self.b.sites() {
            for lam_pos in 0..self.lam.sites() {
                for q_pos in 0..self.q.sites() {
                    let f0 = self.g0[b_pos][lam_pos][q_pos];
                    let f1 = self.g1[b_pos][lam_pos][q_pos];
                    let f2 = self.g2[b_pos][lam_pos][q_pos];

                    // Difference from order variation
                    let d_l = f64::max((f0-f1).abs(), (f0-f2).abs());
                    // Variation of resummation parameters
                    let var_b   = self.b_var(b_pos, lam_pos, q_pos);
                    let var_lam = self.lam_var(b_pos, lam_pos, q_pos);
                    let var_q   = self.q_var(b_pos, lam_pos, q_pos);
                    // Total error
                    let err = d_l + var_b + var_lam + var_q;

                    err_db.push( (err,
                                  b_pos, lam_pos, q_pos, f0,
                                  d_l, var_b, var_lam, var_q) );
                }
                pb.inc(self.q.sites() as u64);
            }
        }
        pb.finish_and_clear();

        // Sort vector using error as a key
        err_db.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        // Point with minimal error
        let e_bar = err_db[0].0;
        // Keep points with E < 3*Ebar
        err_db.retain(|&x| x.0 < 3f64*e_bar);

        let n_err = err_db.len();

        let mut f_mean: f64 = 0f64;
        for e in &err_db {
            f_mean += e.4/(n_err as f64);
        }

        let mut f_sigma: f64 = 0f64;
        for e in &err_db {
            f_sigma += (e.4 - f_mean)*(e.4 - f_mean);
        }
        f_sigma = (f_sigma/(n_err as f64)).sqrt();


        println!("\n{:^40}\n", style("MINIMUM ERROR PARAMETERS").bold().dim());

        println!("{:>10} = {:>20.6}", "Ebar", err_db[0].0);
        println!("{:>10} = {:>20.6}", "d(l)", err_db[0].5);
        println!("{:>10} = {:>20.6}", "V(b)", err_db[0].6);
        println!("{:>10} = {:>20.6}", "V(lam)", err_db[0].7);
        println!("{:>10} = {:>20.6}", "Var(q)", err_db[0].8);
        println!("{:>10} = {:>20.6}", "f", err_db[0].4);
        println!();
        println!("{:>10} = {:>20.6}", "b", self.b.min() + self.b.step()*err_db[0].1 as f64);
        println!("{:>10} = {:>20.6}", "lam", self.lam.min() + self.lam.step()*err_db[0].2 as f64);
        println!("{:>10} = {:>20.6}", "q", self.q.min() + self.q.step()*err_db[0].3 as f64);
        
        // result
        (err_db[0].4, f_mean, f_sigma)
    }
}


fn header() {
    println!("{}", style("                                           ").bold().dim());
    println!("{}", style("   ▄▄▄▄    ▒█████   ██▀███   ██▓  ██████   ").bold().dim());
    println!("{}", style("  ▓█████▄ ▒██▒  ██▒▓██ ▒ ██▒▓██▒▒██    ▒   ").bold().dim());
    println!("{}", style("  ▒██▒ ▄██▒██░  ██▒▓██ ░▄█ ▒▒██▒░ ▓██▄     ").bold().dim());
    println!("{}", style("  ▒██░█▀  ▒██   ██░▒██▀▀█▄  ░██░  ▒   ██▒  ").bold().dim());
    println!("{}", style("  ░▓█  ▀█▓░ ████▓▒░░██▓ ▒██▒░██░▒██████▒▒  ").bold().dim());
    println!("{}", style("  ░▒▓███▀▒░ ▒░▒░▒░ ░ ▒▓ ░▒▓░░▓  ▒ ▒▓▒ ▒ ░  ").bold().dim());
    println!("{}", style("  ▒░▒   ░   ░ ▒ ▒░   ░▒ ░ ▒░ ▒ ░░ ░▒  ░ ░  ").bold().dim());
    println!("{}", style("  ░    ░ ░ ░ ░ ▒    ░░   ░  ▒ ░░  ░  ░     ").bold().dim());
    println!("{}", style("  ░          ░ ░     ░      ░        ░     ").bold().dim());
    println!("{}", style("        ░                                  ").bold().dim());
    println!("{}", style("                                           ").bold().dim());

    println!("{}", style(" Borel resummation with conformal mapping. ").bold().dim());
    println!("{}", style("                                           ").bold().dim());
    println!("{}", style(" Andrey Pikelner  <pikelner@theor.jinr.ru> ").bold().dim());
    println!("{}", style("                                           ").bold().dim());
    println!("{}", style("                                           ").bold().dim());
}


fn main() {


    header();

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

    // 
    // Input from config file only
    // 
    let info_str = match cfg_data.info {
        Some(p) => p,
        None => "".to_string(),
    };
    // input series f = [...]
    let mut f: Vec<f64> = cfg_data.f.iter().map(|&val| val as f64).collect();
    // input ep
    let ep: f64 = cfg_data.ep as f64;
    // input a
    let a: f64 = cfg_data.a as f64;
    let b   = Par3 {l: cfg_data.b[0] as f64, u: cfg_data.b[1] as f64, d: cfg_data.b[2] as usize,};
    let lam = Par3 {l: cfg_data.lam[0] as f64, u: cfg_data.lam[1] as f64, d: cfg_data.lam[2] as usize,};
    let q   = Par3 {l: cfg_data.q[0] as f64, u: cfg_data.q[1] as f64, d: cfg_data.q[2] as usize,};

    let delta_b   = cfg_data.deltab;
    let delta_lam = cfg_data.deltalam;
    let delta_q   = cfg_data.deltaq;



    // workers
    let num_workers = match args.workers {
        Some(n) => {
            rayon::ThreadPoolBuilder::new().num_threads(n).build_global().unwrap();
            n
        },
        None    => rayon::current_num_threads(),
    };
    
    println!("{}   {}", style("[JOBS]").bold().dim(), num_workers);
    println!("{}   {}", style("[INFO]").bold().dim(), &info_str);

    println!("\n{:^40}", style("INPUT PARAMETERS").bold().dim());
    println!("{:>10} = {:<20}", "ep", &ep);
    println!("{:>10} = {:<20}", "a", &a);

    println!("\n{:^40}", style("GRID PARAMETERS").bold().dim());
    println!("{:>10} = [{:3.2}, {:6.2}]/{:<6}", "b", b.min(), b.max(), b.div());
    println!("{:>10} = [{:3.2}, {:6.2}]/{:<6}", "lam", lam.min(), lam.max(), lam.div());
    println!("{:>10} = [{:3.2}, {:6.2}]/{:<6}", "q", q.min(), q.max(), q.div());

    println!("\n{:^40}", style("SCAN PARAMETERS").bold().dim());
    println!("{:>10} = {:<20}", "Delta(b)", &delta_b);
    println!("{:>10} = {:<20}", "Delta(lam)", &delta_lam);
    println!("{:>10} = {:<20}", "Delta(q)", &delta_q);

    println!("\n{:^40}", style("INPUT SERIES").bold().dim());
    for (i, t) in f.iter().enumerate() {
        println!("      x^{:<2} = {:>10.5}", i, t);
    }

    // Set order of integration rules
    let integ_ord = cfg_data.intord.unwrap_or(20);

    let niter = b.sites() * lam.sites() * q.sites();

    println!("\nResummation grid: {}*{}*{} = {}\n", b.sites(), lam.sites(), q.sites(), &niter);


    let mut grid_0 = vec![vec![vec![0f64; q.sites()]; lam.sites()]; b.sites()];
    let mut grid_1 = vec![vec![vec![0f64; q.sites()]; lam.sites()]; b.sites()];
    let mut grid_2 = vec![vec![vec![0f64; q.sites()]; lam.sites()]; b.sites()];

    let b_step   = b.step();
    let lam_step = lam.step();
    let q_step   = q.step();

    // L 
    println!("{} Series order l = {}", style("[1/3]").bold().dim(), f.len() - 1);

    let pb = ProgressBar::new(niter as u64);
    pb.set_style(ProgressStyle::default_bar().template("{spinner:.gray} {bar:40}").unwrap());
    for ib in 0..b.sites() {
        let rsmtr  = boris::resummator::Resummator::new
            (integ_ord as usize,
             a,       // a
             b.min() + b_step*ib as f64,
             ep,      // ep
             f.as_slice()
            );

        for ilam in 0..lam.sites() {
            grid_0[ib][ilam] =
                (0..q.sites()).into_par_iter().map(|iq| rsmtr.point(lam.min() + lam_step*ilam as f64, q.min() + q_step*iq as f64)).collect();

            pb.inc(q.sites() as u64);
        }
    }
    pb.finish_and_clear();

    // Drop last series term
    f.pop();
    // L-1
    println!("{} Series order (l-1) = {}", style("[2/3]").bold().dim(), f.len() - 1);

    let pb = ProgressBar::new(niter as u64);
    pb.set_style(ProgressStyle::default_bar().template("{spinner:.gray} {bar:40}").unwrap());
    for ib in 0..b.sites() {
        let rsmtr  = boris::resummator::Resummator::new
            (integ_ord as usize,
             a,       // a
             b.min() + b_step*ib as f64,
             ep,      // ep
             f.as_slice()
            );

        for ilam in 0..lam.sites() {
            grid_1[ib][ilam] =
                (0..q.sites()).into_par_iter().map(|iq| rsmtr.point(lam.min() + lam_step*ilam as f64, q.min() + q_step*iq as f64)).collect();
            
            pb.inc(q.sites() as u64);
        }

    }
    pb.finish_and_clear();

    // Drop last series term
    f.pop();
    // L-2
    println!("{} Series order (l-2) = {}", style("[3/3]").bold().dim(), f.len() - 1);

    let pb = ProgressBar::new(niter as u64);
    pb.set_style(ProgressStyle::default_bar().template("{spinner:.gray} {bar:40}").unwrap());
    for ib in 0..b.sites() {
        let rsmtr  = boris::resummator::Resummator::new
            (integ_ord as usize,
             a,       // a
             b.min() + b_step*ib as f64,
             ep,      // ep
             f.as_slice()
            );

        for ilam in 0..lam.sites() {
            grid_2[ib][ilam] =
                (0..q.sites()).into_par_iter().map(|iq| rsmtr.point(lam.min() + lam_step*ilam as f64, q.min() + q_step*iq as f64)).collect();

            pb.inc(q.sites() as u64);
        }

    }
    pb.finish_and_clear();


    // We dump grids only if output file name is provided
    if cfg_data.dump.is_some() || args.dump.is_some() {
        let f_dump_state = match args.dump {
            Some(dmp) => File::create(dmp),
            None      => File::create(cfg_data.dump.unwrap()),
        };

        if let Ok(mut f_dump) = f_dump_state {
            
            for ib in 0..(b.div()+1) {
                for ilam in 0..(lam.div()+1) {
                    for iq in 0..(q.div()+1) {
                        writeln!(f_dump, "{}\t{}\t{}\t{}\t{}\t{}",
                                 b.min() + b_step*ib as f64,
                                 lam.min() + lam_step*ilam as f64,
                                 q.min() + q_step*iq as f64,
                                 grid_0[ib][ilam][iq],
                                 grid_1[ib][ilam][iq],
                                 grid_2[ib][ilam][iq],
                        ).unwrap();
                    }
                }
            }
        } else {
            panic!("Error opening DUMP file")
        }
    }

    // Analyzing grids

    let mut err_est = ErrorEstimator::new(grid_0, grid_1, grid_2, b, lam, q);
    let (res_min, mean, sigma) = err_est.scan(delta_b, delta_lam, delta_q);


    println!("\n{}", style("===========================================").bold().dim());
    println!("\n{:^40}\n", style("RESULTS").bold().dim());

    println!("{:>10} = {:>20.6}", style("[F]").bold().dim(), res_min);
    println!("{:>10} = {:>20.6}", style("<F>").bold().dim(), mean);
    println!("{:>10} = {:>20.6}", style("sigma").bold().dim(), sigma);


    // Save results to file
    match cfg_data.out {
        Some(o) => {
            let f_out_state = File::create(o);
            if let Ok(mut f_out) = f_out_state {
                if writeln!(f_out, "{{ {}, {}, {} }}", res_min, mean, sigma).is_ok() {
                    println!("\n{:^40}\n", style("RESULTS SAVED").bold().dim());
                }
            }
        },
        None  => {},
    };
    
}
