#include "Integrator.h"
#include "Interface.h"

#include <sstream>
#include <Rcpp.h>

void cleanup() {
    SimInterface::freeModelLibrary();
}

// Rprintf returns void instead of the std::printf/mex::printf with int -> wrap and return int
int wrapRprint(const char *c,...) {
    va_list(ap);
    va_start(ap, c);
    Rvprintf(c, ap);
    va_end(ap);
    return 0;
}

// [[Rcpp::init]]
void initInterfaceR(DllInfo *dl) {
    SimInterface::setPrintfPtr(&wrapRprint);
    //atExit(cleanup);
}

// [[Rcpp::export(name=".CompiledSimulation_load",rng = false)]]
Rcpp::XPtr<SimInterface> CompiledSimulation_load(const std::string &name) {
    Rcpp::XPtr<SimInterface> sim(new SimInterface(name, SimInterface::enumJacTypes::DENSE_CPP), true);
    if(sim->getErrorFlag())
        throw(Rcpp::exception("Error initializing simulation"));
    else
        return sim;
}

// [[Rcpp::export(name=".CompiledSimulation_run",rng = false)]]
Rcpp::NumericVector CompiledSimulation_run(Rcpp::XPtr<SimInterface> sim, Rcpp::NumericVector p) {
    if(p.size()!=sim->getDimP_free()) {
        std::ostringstream oss;
        oss << "Error parameter vector size mismatch. Simluation has " << sim->getDimP_free() << " free parameters. Got " << p.size() << ".";
        throw(Rcpp::exception(oss.str().c_str()));
    }
    for(unsigned int i=0; i<p.size(); i++)
        sim->wsP[i] = p[i];

    std::vector<double> obs(sim->getDimO()*sim->getDimT());
    bool isSuccess = sim->Simulate(obs.data(),nullptr,nullptr);
    if(isSuccess)
        return Rcpp::wrap(obs);
    else
        throw(Rcpp::exception("Integration failed."));
}

// [[Rcpp::export(name=".CompiledSimulation_getP",rng = false)]]
Rcpp::NumericVector CompiledSimulation_getP(Rcpp::XPtr<SimInterface> sim) {
    Rcpp::NumericVector p(sim->getDimP_free());
    for(unsigned int i=0; i<p.size(); i++)
        p[i] = sim->wsP[i];
    return p;
}

// [[Rcpp::export(name=".CompiledSimulation_getT",rng = false)]]
Rcpp::NumericVector CompiledSimulation_getT(Rcpp::XPtr<SimInterface> sim) {
    Rcpp::NumericVector t(sim->getDimT());
    for(unsigned int i=0; i<t.size(); i++)
        t[i] = sim->wsT[i];
    return t;
}

// [[Rcpp::export(name=".CompiledSimulation_setT",rng = false)]]
void CompiledSimulation_setT(Rcpp::XPtr<SimInterface> sim, Rcpp::NumericVector t) {
    double t_dbl[t.size()];
    for(unsigned int i=0; i<t.size(); i++)
        t_dbl[i] = t[i];
    sim->setT(t.size(), t_dbl);
}
