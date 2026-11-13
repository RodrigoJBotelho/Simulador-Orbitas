// orbitas.cpp
// Problema de dois corpos (satélite em torno da Terra) integrado com quatro
// métodos numéricos, comparados com a solução exata de Kepler.
//
//   d2r/dt2 = -GM r / |r|^3          (plano orbital, 2D)
//
// Métodos:
//   euler       Euler explícito          (ordem 1, a energia CRESCE: a órbita "abre")
//   simpletico  Euler simplético         (ordem 1, energia limitada)
//   verlet      Velocity Verlet          (ordem 2, simplético)
//   rk4         Runge-Kutta clássico     (ordem 4, não simplético)
//
// Condição inicial: periapsis de uma elipse de semieixo a e excentricidade e.
//
// Uso:  ./orbitas [e=0.5] [n_orbitas=5] [passos_por_orbita=200] [--verificar]
// Saídas (no diretório atual):
//   orbita_trajetoria.csv   x_*, y_* em km (para cada método e para a solução exata)
//   orbita_energia.csv      erro relativo de energia x tempo (em órbitas)
//   orbita_convergencia.csv erro de posição após 1 órbita x passos por órbita
//   orbita_convergencia_energia.csv  erro máx. de energia em 1 órbita x passos por órbita
// Com --verificar (0.1 <= e <= 0.8) o programa confere ordens de convergência e conservação de
// energia e retorna código 1 se algo estiver fora do esperado (usado no CI).

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

static const double GM = 3.986004418e14;  // Terra [m^3/s^2]
static const double PI = 3.14159265358979323846;
static const double A_SEMIEIXO = 2.4e7;   // semieixo maior [m] (24 000 km)

struct Estado { double x, y, vx, vy; };

static void acel(const Estado& s, double& ax, double& ay) {
    double r2 = s.x * s.x + s.y * s.y;
    double r3 = r2 * std::sqrt(r2);
    ax = -GM * s.x / r3;
    ay = -GM * s.y / r3;
}

static Estado passo_euler(const Estado& s, double dt) {
    double ax, ay;
    acel(s, ax, ay);
    Estado n = {s.x + dt * s.vx, s.y + dt * s.vy, s.vx + dt * ax, s.vy + dt * ay};
    return n;
}

static Estado passo_simpletico(const Estado& s, double dt) {
    double ax, ay;
    acel(s, ax, ay);
    Estado n;
    n.vx = s.vx + dt * ax;
    n.vy = s.vy + dt * ay;
    n.x = s.x + dt * n.vx;
    n.y = s.y + dt * n.vy;
    return n;
}

static Estado passo_verlet(const Estado& s, double dt) {
    double ax, ay, bx, by;
    acel(s, ax, ay);
    Estado n;
    n.x = s.x + dt * s.vx + 0.5 * dt * dt * ax;
    n.y = s.y + dt * s.vy + 0.5 * dt * dt * ay;
    acel(n, bx, by);
    n.vx = s.vx + 0.5 * dt * (ax + bx);
    n.vy = s.vy + 0.5 * dt * (ay + by);
    return n;
}

static void deriv(const Estado& s, Estado& d) {
    d.x = s.vx;
    d.y = s.vy;
    acel(s, d.vx, d.vy);
}

static Estado soma(const Estado& s, const Estado& d, double h) {
    Estado r = {s.x + h * d.x, s.y + h * d.y, s.vx + h * d.vx, s.vy + h * d.vy};
    return r;
}

static Estado passo_rk4(const Estado& s, double dt) {
    Estado k1, k2, k3, k4;
    deriv(s, k1);
    deriv(soma(s, k1, 0.5 * dt), k2);
    deriv(soma(s, k2, 0.5 * dt), k3);
    deriv(soma(s, k3, dt), k4);
    Estado n;
    n.x = s.x + dt / 6.0 * (k1.x + 2 * k2.x + 2 * k3.x + k4.x);
    n.y = s.y + dt / 6.0 * (k1.y + 2 * k2.y + 2 * k3.y + k4.y);
    n.vx = s.vx + dt / 6.0 * (k1.vx + 2 * k2.vx + 2 * k3.vx + k4.vx);
    n.vy = s.vy + dt / 6.0 * (k1.vy + 2 * k2.vy + 2 * k3.vy + k4.vy);
    return n;
}

typedef Estado (*Passo)(const Estado&, double);

struct Metodo {
    const char* nome;
    Passo passo;
};

static const Metodo METODOS[4] = {
    {"euler", passo_euler},
    {"simpletico", passo_simpletico},
    {"verlet", passo_verlet},
    {"rk4", passo_rk4},
};

// Solução exata de Kepler (elipse percorrida no sentido anti-horário, periapsis em t=0)
struct Orbita { double a, e, n; };

static Estado exato(const Orbita& o, double t) {
    double M = std::fmod(o.n * t, 2.0 * PI);
    double E = (o.e > 0.8) ? PI : M;
    for (int i = 0; i < 60; ++i) {
        double f = E - o.e * std::sin(E) - M;
        double fp = 1.0 - o.e * std::cos(E);
        double d = f / fp;
        E -= d;
        if (std::fabs(d) < 1e-15) break;
    }
    double b = o.a * std::sqrt(1.0 - o.e * o.e);
    double Ed = o.n / (1.0 - o.e * std::cos(E));
    Estado s = {o.a * (std::cos(E) - o.e), b * std::sin(E),
                -o.a * std::sin(E) * Ed, b * std::cos(E) * Ed};
    return s;
}

static double energia(const Estado& s) {
    return 0.5 * (s.vx * s.vx + s.vy * s.vy) - GM / std::sqrt(s.x * s.x + s.y * s.y);
}
static double mom_angular(const Estado& s) { return s.x * s.vy - s.y * s.vx; }

static double erro_pos_km(const Estado& a, const Estado& b) {
    return std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y)) / 1000.0;
}

// erro de posição (km) após 1 órbita com 'passos' passos
static double erro_uma_orbita(const Orbita& o, Passo passo, int passos, double T) {
    Estado s = exato(o, 0.0);
    double dt = T / passos;
    for (int i = 0; i < passos; ++i) s = passo(s, dt);
    return erro_pos_km(s, exato(o, T));
}

// erro relativo MÁXIMO de energia ao longo de 1 órbita com 'passos' passos
static double erro_energia_orbita(const Orbita& o, Passo passo, int passos, double T) {
    Estado s = exato(o, 0.0);
    double E0 = energia(s), mx = 0.0, dt = T / passos;
    for (int i = 0; i < passos; ++i) {
        s = passo(s, dt);
        mx = std::max(mx, std::fabs((energia(s) - E0) / E0));
    }
    return mx;
}

int main(int argc, char** argv) {
    bool verificar = false;
    std::vector<std::string> pos;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--verificar") == 0) verificar = true;
        else pos.push_back(argv[i]);
    }
    double e = pos.size() > 0 ? std::atof(pos[0].c_str()) : 0.5;
    int n_orb = pos.size() > 1 ? std::atoi(pos[1].c_str()) : 5;
    int passos = pos.size() > 2 ? std::atoi(pos[2].c_str()) : 200;
    if (e < 0.0 || e >= 0.95 || n_orb < 1 || passos < 10) {
        std::cerr << "Uso: ./orbitas [0<=e<0.95] [n_orbitas>=1] [passos_por_orbita>=10] [--verificar]\n";
        return 1;
    }

    if (verificar && (e < 0.1 || e > 0.8)) {
        std::cerr << "--verificar so e definido para 0.1 <= e <= 0.8 "
                     "(orbita circular e e muito excentrica tem comportamento de ordem atipico).\n";
        return 2;
    }

    Orbita o;
    o.a = A_SEMIEIXO;
    o.e = e;
    o.n = std::sqrt(GM / (o.a * o.a * o.a));
    double T = 2.0 * PI / o.n;
    double dt = T / passos;
    int total = n_orb * passos;
    Estado s0 = exato(o, 0.0);
    double E0 = energia(s0), L0 = mom_angular(s0);

    std::printf("a=%.0f km  e=%.3f  periodo=%.1f s (%.2f h)  perigeu=%.0f km  apogeu=%.0f km\n",
                o.a / 1e3, e, T, T / 3600.0, o.a * (1 - e) / 1e3, o.a * (1 + e) / 1e3);
    std::printf("passos por orbita=%d  orbitas=%d  dt=%.2f s\n\n", passos, n_orb, dt);

    // ---------- simulação principal ----------
    Estado s[4];
    double maxE[4] = {0, 0, 0, 0}, maxL[4] = {0, 0, 0, 0};
    for (int m = 0; m < 4; ++m) s[m] = s0;

    std::ofstream ft("orbita_trajetoria.csv"), fe("orbita_energia.csv");
    ft << std::setprecision(10);
    fe << std::setprecision(10);
    ft << "x_exato,y_exato";
    fe << "t_orbitas";
    for (int m = 0; m < 4; ++m) {
        ft << ",x_" << METODOS[m].nome << ",y_" << METODOS[m].nome;
        fe << "," << METODOS[m].nome;
    }
    ft << "\n";
    fe << "\n";
    ft << s0.x / 1e3 << "," << s0.y / 1e3;
    for (int m = 0; m < 4; ++m) ft << "," << s0.x / 1e3 << "," << s0.y / 1e3;
    ft << "\n";

    for (int i = 1; i <= total; ++i) {
        double t = i * dt;
        Estado ex = exato(o, t);
        ft << ex.x / 1e3 << "," << ex.y / 1e3;
        fe << (double)i / passos;
        for (int m = 0; m < 4; ++m) {
            s[m] = METODOS[m].passo(s[m], dt);
            double erE = std::fabs((energia(s[m]) - E0) / E0);
            double erL = std::fabs((mom_angular(s[m]) - L0) / L0);
            maxE[m] = std::max(maxE[m], erE);
            maxL[m] = std::max(maxL[m], erL);
            ft << "," << s[m].x / 1e3 << "," << s[m].y / 1e3;
            fe << "," << erE;
        }
        ft << "\n";
        fe << "\n";
    }

    std::printf("%-12s %16s %16s %18s\n", "metodo", "max erro energia", "max erro mom.ang.", "erro posicao final");
    double tf = total * dt;
    for (int m = 0; m < 4; ++m)
        std::printf("%-12s %16.3e %16.3e %15.3f km\n", METODOS[m].nome, maxE[m], maxL[m],
                    erro_pos_km(s[m], exato(o, tf)));

    // ---------- convergência ----------
    const int NS[] = {25, 50, 100, 200, 400, 800, 1600, 3200};
    const int NN = 8;
    double err[4][8], erE[4][8];
    std::ofstream fc("orbita_convergencia.csv"), fk("orbita_convergencia_energia.csv");
    fc << std::setprecision(10);
    fk << std::setprecision(10);
    fc << "passos_por_orbita";
    fk << "passos_por_orbita";
    for (int m = 0; m < 4; ++m) {
        fc << "," << METODOS[m].nome;
        fk << "," << METODOS[m].nome;
    }
    fc << "\n";
    fk << "\n";
    for (int k = 0; k < NN; ++k) {
        fc << NS[k];
        fk << NS[k];
        for (int m = 0; m < 4; ++m) {
            err[m][k] = erro_uma_orbita(o, METODOS[m].passo, NS[k], T);
            erE[m][k] = erro_energia_orbita(o, METODOS[m].passo, NS[k], T);
            fc << "," << err[m][k];
            fk << "," << erE[m][k];
        }
        fc << "\n";
        fk << "\n";
    }
    // ordem observada no erro de energia entre 800 e 1600 passos por órbita
    double ordem[4];
    std::printf("\nOrdem de convergencia observada (erro de energia, 800 -> 1600 passos/orbita):\n");
    for (int m = 0; m < 4; ++m) {
        ordem[m] = std::log(erE[m][5] / erE[m][6]) / std::log(2.0);
        std::printf("  %-12s %.2f\n", METODOS[m].nome, ordem[m]);
    }
    std::printf("\nArquivos gerados: orbita_trajetoria.csv, orbita_energia.csv,\n"
                "                  orbita_convergencia.csv, orbita_convergencia_energia.csv\n");

    if (verificar) {
        const double lo[4] = {0.6, 0.85, 1.85, 3.8};
        const double hi[4] = {1.25, 1.25, 2.15, 4.5};
        bool ok = true;
        for (int m = 0; m < 4; ++m)
            if (!(ordem[m] >= lo[m] && ordem[m] <= hi[m])) {
                std::printf("FALHA: ordem de %s = %.2f fora de [%.2f, %.2f]\n", METODOS[m].nome, ordem[m], lo[m], hi[m]);
                ok = false;
            }
        // a solução exata deve conservar a energia
        double dEx = std::fabs((energia(exato(o, 0.37 * T)) - E0) / E0);
        if (dEx > 1e-9) { std::printf("FALHA: solucao exata nao conserva energia (%.2e)\n", dEx); ok = false; }
        // Euler explícito deve ser pior que Verlet na energia
        if (!(maxE[0] > 10.0 * maxE[2])) { std::printf("FALHA: Euler nao pior que Verlet na energia\n"); ok = false; }
        std::printf(ok ? "VERIFICACAO OK\n" : "VERIFICACAO FALHOU\n");
        return ok ? 0 : 1;
    }
    return 0;
}
