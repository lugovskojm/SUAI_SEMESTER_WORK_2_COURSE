#include "PolyFactorService.h"
#include <sstream>
#include <algorithm>
#include <set>

// ─── Арифметика GF(p) ────────────────────────────────────────────────────────

long long PolyFactorService::mod(long long a) const {
    return ((a % p_) + p_) % p_;
}

long long PolyFactorService::modpow(long long base, long long exp, long long m) const {
    long long result = 1;
    base %= m;
    if (base < 0) base += m;
    while (exp > 0) {
        if (exp & 1) result = result * base % m;
        base = base * base % m;
        exp >>= 1;
    }
    return result;
}

long long PolyFactorService::modinv(long long a) const {
    return modpow(mod(a), p_ - 2, p_);
}

// ─── Многочленная арифметика ──────────────────────────────────────────────────

int PolyFactorService::degree(const Poly& a) const {
    int d = (int)a.size() - 1;
    while (d > 0 && a[d] == 0) d--;
    return d;
}

Poly PolyFactorService::normalize(Poly a) const {
    while (a.size() > 1 && a.back() == 0) a.pop_back();
    if (a.empty()) a.push_back(0);
    long long lc = a.back();
    if (lc != 0 && lc != 1) {
        long long inv = modinv(lc);
        for (auto& c : a) c = mod(c * inv);
    }
    return a;
}

Poly PolyFactorService::polyAdd(const Poly& a, const Poly& b) const {
    Poly res(std::max(a.size(), b.size()), 0);
    for (size_t i = 0; i < a.size(); i++) res[i] = mod(res[i] + a[i]);
    for (size_t i = 0; i < b.size(); i++) res[i] = mod(res[i] + b[i]);
    while (res.size() > 1 && res.back() == 0) res.pop_back();
    return res;
}

Poly PolyFactorService::polySub(const Poly& a, const Poly& b) const {
    Poly res(std::max(a.size(), b.size()), 0);
    for (size_t i = 0; i < a.size(); i++) res[i] = mod(res[i] + a[i]);
    for (size_t i = 0; i < b.size(); i++) res[i] = mod(res[i] - b[i] + p_);
    while (res.size() > 1 && res.back() == 0) res.pop_back();
    return res;
}

Poly PolyFactorService::polyMul(const Poly& a, const Poly& b) const {
    if (a.empty() || b.empty()) return {0};
    Poly res(a.size() + b.size() - 1, 0);
    for (size_t i = 0; i < a.size(); i++)
        for (size_t j = 0; j < b.size(); j++)
            res[i + j] = mod(res[i + j] + a[i] * b[j]);
    while (res.size() > 1 && res.back() == 0) res.pop_back();
    return res;
}

Poly PolyFactorService::polyScale(const Poly& a, long long s) const {
    Poly res = a;
    for (auto& c : res) c = mod(c * s);
    return res;
}

// Возвращает остаток от деления a на b
Poly PolyFactorService::polyMod(const Poly& a, const Poly& b) const {
    Poly rem = a;
    while (rem.size() > 1 && rem.back() == 0) rem.pop_back();
    int db = degree(b);
    long long lc_inv = modinv(b[db]);
    while (degree(rem) >= db && !(rem.size() == 1 && rem[0] == 0)) {
        int dr = degree(rem);
        if (dr < db) break;
        long long coef = mod(rem[dr] * lc_inv);
        int shift = dr - db;
        for (int i = 0; i <= db; i++) {
            int idx = i + shift;
            if (idx < (int)rem.size())
                rem[idx] = mod(rem[idx] - b[i] * coef + p_);
        }
        while (rem.size() > 1 && rem.back() == 0) rem.pop_back();
    }
    return rem;
}

// Возвращает частное от деления a на b
Poly PolyFactorService::polyDiv(const Poly& a, const Poly& b) const {
    Poly rem = a;
    while (rem.size() > 1 && rem.back() == 0) rem.pop_back();
    int db = degree(b);
    long long lc_inv = modinv(b[db]);
    Poly quot;
    while (degree(rem) >= db && !(rem.size() == 1 && rem[0] == 0)) {
        int dr = degree(rem);
        if (dr < db) break;
        long long coef = mod(rem[dr] * lc_inv);
        int shift = dr - db;
        if ((int)quot.size() <= shift) quot.resize(shift + 1, 0);
        quot[shift] = mod(quot[shift] + coef);
        for (int i = 0; i <= db; i++) {
            int idx = i + shift;
            if (idx < (int)rem.size())
                rem[idx] = mod(rem[idx] - b[i] * coef + p_);
        }
        while (rem.size() > 1 && rem.back() == 0) rem.pop_back();
    }
    if (quot.empty()) quot.push_back(0);
    while (quot.size() > 1 && quot.back() == 0) quot.pop_back();
    return quot;
}

Poly PolyFactorService::polyGCD(Poly a, Poly b) const {
    a = normalize(a);
    b = normalize(b);
    while (!(b.size() == 1 && b[0] == 0)) {
        Poly r = polyMod(a, b);
        a = b;
        b = normalize(r);
    }
    return normalize(a);
}

Poly PolyFactorService::polyDeriv(const Poly& f) const {
    if (f.size() <= 1) return {0};
    Poly d(f.size() - 1);
    for (size_t i = 1; i < f.size(); i++)
        d[i - 1] = mod(f[i] * (long long)i);
    while (d.size() > 1 && d.back() == 0) d.pop_back();
    return d;
}

Poly PolyFactorService::polyPowMod(const Poly& base, long long exp,
                                    const Poly& modulus) const {
    Poly result = {1};
    Poly b = polyMod(base, modulus);
    while (exp > 0) {
        if (exp & 1) result = polyMod(polyMul(result, b), modulus);
        b = polyMod(polyMul(b, b), modulus);
        exp >>= 1;
    }
    return result;
}

// ─── Square-free decomposition ────────────────────────────────────────────────

std::vector<std::pair<Poly, int>>
PolyFactorService::squareFreeDecompose(const Poly& f) const {
    std::vector<std::pair<Poly, int>> result;

    Poly fd = polyDeriv(f);

    // Если производная = 0, то f = g^p (характеристика поля)
    if (fd.size() == 1 && fd[0] == 0) {
        // Понижаем степень: f(x) = g(x^p) => извлекаем g
        // Коэффициенты с индексами кратными p
        Poly g;
        int deg = degree(f);
        for (int i = 0; i <= deg; i += (int)p_) {
            g.push_back(f[i]);
        }
        while (g.size() > 1 && g.back() == 0) g.pop_back();
        // Рекурсивно раскладываем g, умножая кратности на p
        auto sub = squareFreeDecompose(g);
        for (auto& [poly, mult] : sub)
            result.push_back({poly, mult * (int)p_});
        return result;
    }

    Poly g = polyGCD(f, fd);

    if (degree(g) == 0) {
        // f уже бесквадратное
        result.push_back({normalize(f), 1});
        return result;
    }

    Poly w = polyDiv(f, g);  // w = f / gcd(f, f')
    int i = 1;

    while (degree(w) > 0) {
        Poly y = polyGCD(w, g);
        Poly factor = polyDiv(w, y);  // factor = w / y
        if (degree(factor) > 0)
            result.push_back({normalize(factor), i});
        w = y;
        g = polyDiv(g, y);
        i++;
    }

    // Оставшаяся часть g имеет все множители с кратностями кратными p
    if (degree(g) > 0) {
        Poly gp;
        int deg = degree(g);
        for (int j = 0; j <= deg; j += (int)p_) {
            gp.push_back(g[j]);
        }
        while (gp.size() > 1 && gp.back() == 0) gp.pop_back();
        auto sub = squareFreeDecompose(gp);
        for (auto& [poly, mult] : sub)
            result.push_back({poly, mult * (int)p_});
    }

    return result;
}

// ─── Алгоритм Берлекэмпа ─────────────────────────────────────────────────────

std::vector<Poly> PolyFactorService::berlekamp(const Poly& f) const {
    int n = degree(f);
    if (n <= 1) return {normalize(f)};

    // Строим матрицу B: строка i = коэффициенты x^(i*p) mod f
    std::vector<std::vector<long long>> B(n, std::vector<long long>(n, 0));
    Poly xp = {0, 1}; // x
    xp = polyPowMod(xp, p_, f); // x^p mod f

    // Строка 0: x^0 mod f = 1
    Poly cur = {1};
    for (int i = 0; i < n; i++) {
        // cur = x^(i*p) mod f
        Poly row = cur;
        for (int j = 0; j < n && j < (int)row.size(); j++)
            B[i][j] = row[j];
        // следующая строка: умножаем на x^p mod f
        cur = polyMod(polyMul(cur, xp), f);
    }

    // Q = B^T - I
    std::vector<std::vector<long long>> Q(n, std::vector<long long>(n, 0));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            Q[i][j] = (i == j) ? mod(B[j][i] - 1) : B[j][i];

    // Приведение к ступенчатому виду (Гаусс) для нахождения ядра
    std::vector<int> pivotRow(n, -1);
    int rank = 0;
    for (int col = 0; col < n && rank < n; col++) {
        int row = -1;
        for (int i = rank; i < n; i++)
            if (Q[i][col] != 0) { row = i; break; }
        if (row == -1) continue;
        std::swap(Q[rank], Q[row]);
        pivotRow[col] = rank;
        long long inv = modinv(Q[rank][col]);
        for (int j = 0; j < n; j++) Q[rank][j] = mod(Q[rank][j] * inv);
        for (int i = 0; i < n; i++) {
            if (i != rank && Q[i][col] != 0) {
                long long c = Q[i][col];
                for (int j = 0; j < n; j++)
                    Q[i][j] = mod(Q[i][j] - c * Q[rank][j]);
            }
        }
        rank++;
    }

    // Извлекаем базис ядра
    std::vector<Poly> basis;
    std::set<int> pivotCols;
    for (int col = 0; col < n; col++)
        if (pivotRow[col] >= 0) pivotCols.insert(col);

    for (int col = 0; col < n; col++) {
        if (pivotCols.count(col)) continue;
        Poly ker(n, 0);
        ker[col] = 1;
        for (int c2 = 0; c2 < n; c2++) {
            if (!pivotCols.count(c2)) continue;
            int r = pivotRow[c2];
            ker[c2] = mod(-Q[r][col]);
        }
        while (ker.size() > 1 && ker.back() == 0) ker.pop_back();
        basis.push_back(ker);
    }

    // Если ядро тривиально — многочлен неприводим
    if (basis.empty()) return {normalize(f)};

    // Расщепляем f с помощью базисных векторов
    return splitFactors({normalize(f)}, basis);
}

std::vector<Poly> PolyFactorService::splitFactors(
    const std::vector<Poly>& initial,
    const std::vector<Poly>& basis) const
{
    std::vector<Poly> factors = initial;

    for (const Poly& h : basis) {
        std::vector<Poly> newFactors;
        for (const Poly& u : factors) {
            if (degree(u) <= 1) {
                newFactors.push_back(u);
                continue;
            }
            bool wasSplit = false;
            for (long long s = 0; s < p_; s++) {
                // Вычисляем h - s
                Poly hs = h;
                if (hs.empty()) hs.push_back(0);
                while ((int)hs.size() < 1) hs.push_back(0);
                hs[0] = mod(hs[0] - s);
                while (hs.size() > 1 && hs.back() == 0) hs.pop_back();

                if (hs.size() == 1 && hs[0] == 0) continue;

                Poly g = polyGCD(u, hs);
                int dg = degree(g);
                int du = degree(u);

                if (dg > 0 && dg < du) {
                    Poly q = polyDiv(u, g);
                    newFactors.push_back(normalize(g));
                    newFactors.push_back(normalize(q));
                    wasSplit = true;
                    break;
                }
            }
            if (!wasSplit) newFactors.push_back(u);
        }
        factors = newFactors;
    }

    // Нормализуем и убираем единицы
    std::vector<Poly> result;
    for (auto& fac : factors) {
        fac = normalize(fac);
        if (degree(fac) > 0) result.push_back(fac);
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

// ─── Полная факторизация ──────────────────────────────────────────────────────

std::vector<std::pair<Poly, int>>
PolyFactorService::fullFactorize(const Poly& f) const {
    std::vector<std::pair<Poly, int>> result;

    auto sqfParts = squareFreeDecompose(f);

    for (auto& [part, mult] : sqfParts) {
        if (degree(part) <= 0) continue;
        // Применяем Берлекэмп к каждой бесквадратной части
        auto irredFactors = berlekamp(part);
        for (auto& fac : irredFactors) {
            if (degree(fac) > 0)
                result.push_back({normalize(fac), mult});
        }
    }

    // Сортируем по степени
    std::sort(result.begin(), result.end(),
        [](const std::pair<Poly,int>& a, const std::pair<Poly,int>& b) {
            return a.first.size() < b.first.size();
        });

    return result;
}

// ─── Строковое представление ──────────────────────────────────────────────────

std::string PolyFactorService::polyToString(const Poly& poly) const {
    int deg = degree(poly);
    if (deg == 0) return std::to_string(poly[0]);

    std::string s;
    bool first = true;
    for (int i = deg; i >= 0; i--) {
        long long c = (i < (int)poly.size()) ? poly[i] : 0;
        if (c == 0) continue;
        if (!first) s += " + ";
        if (i == 0) {
            s += std::to_string(c);
        } else if (i == 1) {
            if (c != 1) s += std::to_string(c);
            s += "x";
        } else {
            if (c != 1) s += std::to_string(c);
            s += "x^" + std::to_string(i);
        }
        first = false;
    }
    return s.empty() ? "0" : s;
}

// ─── Точка входа ─────────────────────────────────────────────────────────────

std::string PolyFactorService::factorize(long long p, const Poly& poly) {
    p_ = p;
    Poly f = normalize(poly);
    std::ostringstream oss;
    oss << "Многочлен: " << polyToString(f) << " над GF(" << p << ")\n\n";

    if (degree(f) == 0) {
        oss << "Константа, факторизация невозможна.";
        return oss.str();
    }

    auto allFactors = fullFactorize(f);

    bool isIrreducible = (allFactors.size() == 1 &&
                          allFactors[0].second == 1 &&
                          degree(allFactors[0].first) == degree(f));

    if (isIrreducible) {
        oss << "Многочлен неприводим над GF(" << p << ").";
    } else {
        oss << "Разложение на множители:\n";
        oss << polyToString(f) << " = ";
        bool first = true;
        for (auto& [fac, mult] : allFactors) {
            for (int i = 0; i < mult; i++) {
                if (!first) oss << " * ";
                oss << "(" << polyToString(fac) << ")";
                first = false;
            }
        }
        oss << "\n\nМножители:\n";
        int idx = 1;
        for (auto& [fac, mult] : allFactors) {
            oss << "  f" << idx++ << "(x) = " << polyToString(fac);
            if (mult > 1) oss << "  (кратность " << mult << ")";
            oss << "\n";
        }
    }
    return oss.str();
}
