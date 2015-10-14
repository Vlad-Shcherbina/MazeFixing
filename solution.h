#include <assert.h>
#include <vector>
#include <string>
#include <random>
#include <array>
#include <iomanip>
#include <algorithm>
#include <thread>
#include <chrono>
#include <functional>
#include <sstream>

#include "pretty_printing.h"

using namespace std;


#define debug(x) cerr << #x " = " << (x) << endl
#define debug2(x, y) cerr << #x " = " << (x) << ", " #y " = " << y << endl
#define debug3(x, y, z) \
    cerr << #x " = " << (x) << ", " #y " = " << y << ", " #z " = " << z << endl


int W;
int H;
int F;
vector<string> maze;


typedef pair<int, int> Coord;

const int DIRS[4][2] = {{1, 0}, {0, -1}, {-1, 0}, {0, 1}};

typedef unsigned int Dir;

struct Spin {
    int x, y;
    Dir dir;

    Spin(): x(-100), y(-100), dir(-100) {}
    Spin(int x, int y, Dir dir): x(x), y(y), dir(dir % 4) {}

    Coord get_coord() const {
        return {x, y};
    }

    void advance() {
        x += DIRS[dir][0];
        y += DIRS[dir][1];
    }

    void turn_left() {
        dir += 1;
        dir %= 4;
    }
    void turn_right() {
        dir += 3;
        dir %= 4;
    }
};

std::ostream& operator<<(std::ostream &out, Spin s) {
    out << "Spin(" << s.x << ", " << s.y << ", " << s.dir << ")";
    return out;
}

vector<Spin> find_enters() {
    vector<Spin> result;

    for (int y = 1; y < H - 1; y++) {
        for (int x = 1; x < W - 1; x++) {
            if (maze[y][x] != '.') {
                for (Dir dir = 0; dir < 4; dir++) {
                    int dx = DIRS[dir][0];
                    int dy = DIRS[dir][1];
                    if (maze[y - dy][x - dx] == '.') {
                        result.emplace_back(x, y, dir);
                    }
                }
            }
        }
    }
    return result;
}


class SolutionChecker {
private:
    vector<Coord> current_path;
public:
    set<Coord> covered;
    void trace(Spin spin) {
        int original_size = current_path.size();

        // Tail recursion
        while (true) {
            Coord pos = spin.get_coord();

            // PERF: potentially quadratic
            if (find(current_path.begin(), current_path.end(), pos) !=
                current_path.end()) {
                break;
            }

            char c = maze[pos.second][pos.first];

            if (c == '.') {
                // debug(current_path);
                copy(current_path.begin(), current_path.end(),
                     inserter(covered, covered.end()));
                break;
            }

            current_path.push_back(pos);

            switch (c) {
            case 'E':
                for (Dir dir = 0; dir < 4; dir++) {
                    Spin s2 = spin;
                    s2.dir = dir;
                    s2.advance();
                    if (dir < 3)
                        trace(s2);
                    else
                        spin = s2;  // last one is tail call
                }
                continue;
            case 'S':
                spin.advance();
                continue;
            case 'L':
                spin.turn_left();
                spin.advance();
                continue;
            case 'R':
                spin.turn_right();
                spin.advance();
                continue;
            case 'U':
                spin.turn_left();
                spin.turn_left();
                spin.advance();
                continue;
            }
        }

        assert(current_path.size() >= original_size);
        current_path.resize(original_size);
    }
};


class MazeFixing {
public:
    vector<string> improve(vector<string> maze, int F) {
        ::H = maze.size();
        ::W = maze.front().size();
        ::F = F;
        ::maze = maze;

        debug3(W, H, F);

        int total_cells = 0;
        for (const string &row : maze)
            total_cells += row.size() - count(row.begin(), row.end(), '.');

        SolutionChecker checker;
        for (Spin enter : find_enters()) {
            checker.trace(enter);
        }

        for (int i = 0; i < H; i++) {
            for (int j = 0; j < W; j++) {
                if (checker.covered.count({j, i}) > 0)
                    cerr << "* ";
                else
                    cerr << ". ";
            }
            cerr << endl;
        }

        double predicted_score = 1.0 * checker.covered.size() / total_cells;
        debug(predicted_score);

        vector<string> result;
        return result;
    }

    static string format_result(int x, int y, char label) {
        ostringstream out;
        out << y << " " << x << " " << label;
        return out.str();
    }
};
