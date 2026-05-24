#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <limits>
#include <fstream>
#include "PolyFactorService.h"

static void clearScreen() {
    std::cout << "\x1B[2J\x1B[H"; // ANSI clear and move cursor home
}

static long long readLong(const std::string& prompt) {
    while (true) {
        std::cout << prompt;
        std::string line;
        if (!std::getline(std::cin, line)) return 0;
        std::istringstream in(line);
        long long v;
        if (in >> v) return v;
        std::cout << "Неверный ввод, повторите.\n";
    }
}

static std::vector<long long> readCoeffs(const std::string& prompt) {
    std::cout << prompt;
    std::string line;
    std::getline(std::cin, line);
    std::istringstream in(line);
    std::vector<long long> poly;
    long long v;
    while (in >> v) poly.push_back(v);
    return poly;
}

int main() {
    PolyFactorService svc;
    while (true) {
        clearScreen();
        std::cout << "=== PolyFactorGF — TUI ===\n";
        std::cout << "1) Факторизовать многочлен\n";
        std::cout << "2) Загрузить из файла (каждая строка: p coeffs)\n";
        std::cout << "3) Выход\n";
        std::cout << "Выберите опцию: ";
        std::string opt;
        if (!std::getline(std::cin, opt)) break;
        if (opt == "1") {
            long long p = readLong("Введите простое p: ");
            if (p <= 1) { std::cout << "p должно быть >=2. Нажмите Enter...\n"; std::cin.get(); continue; }
            std::cout << "Введите коэффициенты (от младшего к старшему) через пробел: ";
            std::vector<long long> poly = readCoeffs("");
            if (poly.empty()) { std::cout << "Пустой многочлен. Нажмите Enter...\n"; std::cin.get(); continue; }
            std::string out = svc.factorize(p, poly);
            std::cout << "\n" << out << "\n";
            std::cout << "Нажмите Enter для возврата в меню...";
            std::cin.get();
        } else if (opt == "2") {
            std::cout << "Введите путь к файлу: ";
            std::string path;
            if (!std::getline(std::cin, path)) break;
            std::ifstream file(path);
            if (!file) { std::cout << "Не удалось открыть файл. Нажмите Enter...\n"; std::cin.get(); continue; }
            std::string line;
            while (std::getline(file, line)) {
                if (line.empty()) continue;
                std::istringstream in(line);
                long long p;
                if (!(in >> p)) continue;
                std::vector<long long> poly;
                long long v;
                while (in >> v) poly.push_back(v);
                if (poly.empty()) continue;
                std::cout << "\n--- p=" << p << " coeffs=";
                for (size_t i = 0; i < poly.size(); ++i) {
                    if (i) std::cout << ' ';
                    std::cout << poly[i];
                }
                std::cout << " ---\n";
                std::cout << svc.factorize(p, poly) << "\n";
            }
            std::cout << "Готово. Нажмите Enter...";
            std::cin.get();
        } else if (opt == "3") break;
        else { std::cout << "Неверный выбор. Нажмите Enter...\n"; std::cin.get(); }
    }
    return 0;
}
