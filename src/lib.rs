use pyo3::prelude::*;
use rand::Rng;

fn get_mu(a: &Vec<f64>, r2: f64) -> usize {
    let mut total: f64 = 0.0;
    let a0: f64 = a.iter().sum();
    let r2a0: f64 = r2 * a0;
    for (i, ai) in a.iter().enumerate() {
        total += ai;
        if total >= r2a0 {
            return i;
        }
    }
    return 0;
}

/// Simulates a stochastic kinetic process using the Gillespie algorithm.
///
/// Parameters
/// ----------
/// yi : List[int]
///     Initial number of molecules of each species.
/// c : List[float]
///     Reaction rate constants.
/// update_y : List[List[int]]
///     Change in number of molecules of each species for each reaction.
/// func : Callable[[List[int]], List[float]]
///    Function that returns the propensity for each reaction.
/// stop_time : float
///    Time at which to stop the simulation.
///
/// Returns
/// -------
/// Tuple [List[float], List[List[int]]]
///    List of times and list of number of molecules of each species at each time.
///
#[pyfunction]
fn gillespie(
    py: Python,
    yi: Vec<i64>,
    c: Vec<f64>,
    update_y: Vec<Vec<i64>>,
    func: PyObject,
    stop_time: f64,
) -> PyResult<(Vec<f64>, Vec<Vec<i64>>)> {
    let mut times: Vec<f64> = Vec::new();
    let mut species: Vec<Vec<i64>> = Vec::new();

    let mut rng = rand::thread_rng();
    let mut t: f64 = 0.0;
    let mut r1: f64;
    let mut r2: f64;
    let mut a: Vec<f64>;
    let mut a0: f64;
    let mut mu: usize;
    let mut y = yi.clone();

    times.push(t);
    species.push(y.to_vec());
    while t < stop_time {
        let h: Vec<f64> = func
            .call1(py, (y.clone(),))
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyTypeError, _>(e.to_string()))?
            .extract(py)?;
        a = c.iter().zip(h.iter()).map(|(x, y)| x * y).collect();
        a0 = a.iter().sum();
        if a0 == 0.0 {
            break;
        }
        r1 = rng.gen();
        r2 = rng.gen();
        t += (1.0 / r1).ln() / a0;
        mu = get_mu(&a, r2);
        for i in 0..y.len() {
            y[i] += update_y[mu][i];
        }

        times.push(t);
        species.push(y.to_vec());
    }
    Ok((times, species))
}

/// A Python module for stochastic kinetic modeling implemented in Rust.
#[pymodule]
fn gillespy(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(gillespie, m)?)?;
    Ok(())
}
