#include <assert.h>
#include <vector>
#include <string>
#include <random>
#include <array>
#include <algorithm>
#include <thread>
#include <chrono>
#include <functional>
#include <sstream>
#include <sys/time.h>
#include <time.h>

#include "pretty_printing.h"

using namespace std;


#define debug(x) \
    cerr << #x " = " << (x) << endl
#define debug2(x, y) \
    cerr << #x " = " << (x) \
    << ", " #y " = " << (y) << endl
#define debug3(x, y, z) \
    cerr << #x " = " << (x) \
    << ", " #y " = " << (y) \
    << ", " #z " = " << (z) << endl


const double TIME_LIMIT = 9.5;
double get_time() {
  timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec * 1e-6;
}


const int MAX_SIZE = 80;

int W;
int H;
int F;
vector<string> maze, original_maze;

typedef pair<int, int> Coord;


typedef unsigned short PackedCoord;

PackedCoord pack(int x, int y) {
    return x + 100 * y;
}
PackedCoord pack(Coord coord) {
    return pack(coord.first, coord.second);
}
int unpack_x(PackedCoord p) {
    return p % 100;
}
int unpack_y(PackedCoord p) {
    return p / 100;
}
Coord unpack(PackedCoord p) {
    return {unpack_x(p), unpack_y(p)};
}


typedef PackedCoord Vertex;
const Vertex ENTER = 0;
const Vertex EXIT = 9999;
const Vertex INVALID = 42424;


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
            default:
                assert(false);
            }
        }

        assert(current_path.size() >= original_size);
        current_path.resize(original_size);
    }
};


vector<pair<Vertex, Spin>> enumerate_starting_points() {
    vector<pair<Vertex, Spin>> result;
    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            char c = maze[y][x];
            if (c == '.')
                continue;
            for (Dir d = 0; d < 4; d++) {
                int dx = DIRS[d][0];
                int dy = DIRS[d][1];
                char pc = maze[y - dy][x - dx];
                if (pc == '.') {
                    result.emplace_back(ENTER, Spin(x, y, d));
                }
                if (c == 'E') {
                    result.emplace_back(pack(x, y), Spin(x + dx, y + dy, d));
                }
            }
        }
    }
    return result;
}


vector<Spin> find_enters() {
    vector<Spin> result;
    for (auto p : enumerate_starting_points()) {
        if (p.first == ENTER)
            result.push_back(p.second);
    }
    return result;
}


class CellSet {
private:
    mutable vector<PackedCoord> cells;
    mutable bool sorted = true;
public:
    void clear() {
        *this = CellSet();
    }

    void ensure_sorted() const {
        if (!sorted) {
            sort(cells.begin(), cells.end());
            auto p = unique(cells.begin(), cells.end());
            cells.erase(p, cells.end());
            sorted = true;
        }
    }

    int size() const {
        ensure_sorted();
        return cells.size();
    }

    void add_cell(PackedCoord p) {
        cells.push_back(p);
        sorted = false;
    }

    void merge(const CellSet &other) {
        for (PackedCoord p : other.cells)
            add_cell(p);
    }

    bool overlaps(const CellSet &other) const {
        ensure_sorted();
        other.ensure_sorted();
        vector<PackedCoord> intersection;
        set_intersection(
            cells.begin(), cells.end(),
            other.cells.begin(), other.cells.end(),
            back_inserter(intersection));
        return !intersection.empty();
    }

    bool contains(PackedCoord p) const {
        ensure_sorted();
        return binary_search(cells.begin(), cells.end(), p);
    }

    const vector<PackedCoord>& get_cells() const {
        return cells;
    }

    friend ostream& operator<<(ostream &out, const CellSet &cs);
};

ostream& operator<<(ostream &out, const CellSet &cs) {
    for (int y = 0; y < H; y++) {
        out << endl;
        for (int x = 0; x < W; x++) {
            if (cs.contains(pack(x, y)))
                out << "* ";
            else
                out << ". ";
        }
    }
    return out;
}


// TODO: also return covered cells
Vertex trace(Spin s, CellSet &cell_set) {
    vector<PackedCoord> visited;
    while (true) {
        char c = maze[s.y][s.x];
        if (c == '.') {
            return EXIT;
        }
        if (c == 'E') {
            return pack(s.x, s.y);
        }
        PackedCoord pos = pack(s.x, s.y);
        if (find(visited.begin(), visited.end(), pos) != visited.end()) {
            return INVALID;
        }
        visited.push_back(pos);
        cell_set.add_cell(pos);
        switch (c) {
        case 'S':
            s.advance();
            break;
        case 'L':
            s.turn_left();
            s.advance();
            break;
        case 'R':
            s.turn_right();
            s.advance();
            break;
        case 'U':
            s.turn_left();
            s.turn_left();
            s.advance();
            break;
        default:
            assert(false);
        }
    }
}


int eval_maze() {
    SolutionChecker checker;
    for (Spin enter : find_enters()) {
        checker.trace(enter);
    }
    return checker.covered.size();
}


struct Edge {
private:
    Vertex _from;
    Vertex _to;
    map<PackedCoord, char> edits;
    Spin start_spin;
    CellSet _path_cells;

public:
    Edge(
        Vertex from,
        Vertex to,
        const map<PackedCoord, char> &edits,
        Spin start_spin,
        const vector<PackedCoord> path)
        : _from(from), _to(to), edits(edits), start_spin(start_spin)
    {
        for (PackedCoord p : path)
            _path_cells.add_cell(p);
    }

    bool contradicts(const Edge &other) const {
        for (auto kv : edits)
            if (other._path_cells.contains(kv.first))
                return true;
        for (auto kv : other.edits)
            if (_path_cells.contains(kv.first))
                return true;
        return false;
    }

    int from() const { return _from; }
    int to() const { return _to; }
    int path_length() const { return _path_cells.size(); }
    const CellSet& path_cells() const { return _path_cells; }
};


class EdgeFinder {
public:
    vector<Edge> edges;
private:
    vector<PackedCoord> current_path;
    map<PackedCoord, char> edits;
    Spin start_spin;
    Vertex from;

    void add_edge(Vertex to) {
        edges.emplace_back(from, to, edits, start_spin, current_path);
    }

    void trace(Spin s, int num_edits) {
        char c = maze[s.y][s.x];
        if (c == '.' || c == 'E') {
            if (num_edits == 0) {
                Vertex to = c == '.' ? EXIT : pack(s.x, s.y);
                add_edge(to);
            }
            return;
        }

        PackedCoord pos = pack(s.x, s.y);
        if (find(current_path.begin(), current_path.end(), pos) !=
            current_path.end()) {
            return;
        }

        current_path.push_back(pos);

        Dir new_dir;
        switch (c) {
        case 'S':
            new_dir = s.dir;
            break;
        case 'L':
            new_dir = s.dir + 1;
            break;
        case 'R':
            new_dir = s.dir + 3;
            break;
        case 'U':
            new_dir = s.dir + 2;
            break;
        default:
            assert(false);
        }
        new_dir %= 4;

        Spin new_spin = s;
        new_spin.dir = new_dir;
        new_spin.advance();
        trace(new_spin, num_edits);

        if (num_edits > 0) {
            for (Dir d = 0; d < 4; d++) {
                if (d != new_dir) {
                    edits[pos] = "SLUR"[(4 + d - s.dir) % 4];
                    new_spin = s;
                    new_spin.dir = d;
                    new_spin.advance();
                    trace(new_spin, num_edits - 1);
                }
            }
            edits.erase(pos);
        }

        assert(current_path.back() == pos);
        current_path.pop_back();
    }

public:
    void trace(Vertex from, Spin s, int num_edits) {
        this->from = from;
        this->start_spin = s;
        trace(s, num_edits);
    }
};


class CellMultiSet {
private:
    vector<int> counts = vector<int>(pack(MAX_SIZE, MAX_SIZE), 0);
    int _size = 0;

public:
    int size() const { return _size; }

    void add(PackedCoord p) {
        if (counts[p] == 0)
            _size++;
        counts[p]++;
    }
    void subtract(PackedCoord p) {
        assert(counts[p] > 0);
        counts[p]--;
        if (counts[p] == 0)
            _size--;
    }

    void add(const CellSet& cs) {
        for (PackedCoord p : cs.get_cells())
            add(p);
    }
    void subtract(const CellSet& cs) {
        for (PackedCoord p : cs.get_cells())
            subtract(p);
    }

    friend ostream& operator<<(ostream& out, const CellMultiSet &cms);
};

ostream& operator<<(ostream& out, const CellMultiSet &cms) {
    for (int i = 0; i < H; i++) {
        out << endl;
        for (int j = 0; j < W; j++) {
            int c = cms.counts[pack(j, i)];
            if (c == 0)
                out << ". ";
            else {
                if (c > 9) c = 0;
                out << c << ' ';
            }
        }
    }
    return out;
}


typedef map<Vertex, vector<Edge>> Graph;

class GraphChecker {
private:
    const Graph &graph;
    vector<Vertex> current_path;
    vector<const Edge*> current_edges;
public:
    CellSet covered;
    GraphChecker(const Graph &graph) : graph(graph) {}

    int get_covered_area() const {
        return covered.size();
    }

    void trace(Vertex start) {
        assert(current_path.empty());
        assert(current_edges.empty());
        current_path.push_back(start);
        rec();
    }

    void rec() {
        Vertex v = current_path.back();
        if (v == EXIT) {
            for (Vertex p : current_path)
                if (p != ENTER && p != EXIT)
                    covered.add_cell(p);
            for (auto e : current_edges)
                covered.merge(e->path_cells());
            return;
        }

        auto kv = graph.find(v);
        if (kv == graph.end()) {
            return;
        }

        for (const auto& e : kv->second) {
            Vertex w = e.to();
            if (find(current_path.begin() + 1, current_path.end(), w) !=
                current_path.end())
                continue;
            bool overlaps = false;
            for (const Edge* pe : current_edges) {
                if (e.path_cells().overlaps(pe->path_cells())) {
                    overlaps = true;
                    break;
                }
            }
            if (overlaps)
                continue;

            current_path.push_back(w);
            current_edges.push_back(&e);

            rec();

            assert(current_path.back() == w);
            current_path.pop_back();
            assert(current_edges.back() == &e);
            current_edges.pop_back();
        }
    }
};

class MazeFixing {
public:
    vector<string> improve(vector<string> maze, int F) {
        ::H = maze.size();
        ::W = maze.front().size();
        ::F = F;
        ::maze = maze;
        ::original_maze = maze;

        double start = get_time();

        debug3(W, H, F);

        int total_cells = 0;
        for (const string &row : maze)
            total_cells += row.size() - count(row.begin(), row.end(), '.');

        for (const string &row : maze) {
            for (char c : row) {
                cerr << c << ' ';
            }
            cerr << endl;
        }

        int num_e = 0;
        for (const string &row : maze)
            num_e += count(row.begin(), row.end(), 'E');
        debug(num_e);

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

        {
            map<pair<bool, bool>, int> cnt;
            for (auto p : enumerate_starting_points()) {
                Vertex start = p.first;
                CellSet cell_set;
                Vertex end = trace(p.second, cell_set);
                if (end != INVALID && start != end) {
                    cnt[{start == ENTER, end == EXIT}]++;
                }
            }
            debug(cnt);
        }

        vector<Edge> candidate_edges;

        Graph graph;

        for (int num_edits = 0; num_edits < 1; num_edits++) {
            debug(num_edits);

            map<pair<int, int>, int> cnt;
            map<int, int> lengths;

            EdgeFinder ef;

            for (auto p : enumerate_starting_points()) {
                Vertex start = p.first;
                ef.trace(start, p.second, num_edits);
            }
            for (Edge &e : ef.edges) {
                if (num_edits == 0)
                    graph[e.from()].push_back(e);

                if (e.from() == ENTER && e.to() == EXIT)
                    candidate_edges.push_back(e);

                cnt[{e.from() == ENTER, e.to() == EXIT}]++;
                lengths[e.path_length()]++;
            }
            debug(cnt);
            debug(cnt.size());
            debug(lengths);
        }
        debug(candidate_edges.size());

        double predicted_score = 1.0 * checker.covered.size() / total_cells;
        debug(predicted_score);

        for (int i = 0; i < 1; i++)
            eval_maze();

        for (int i = 0; i < 1; i++) {
            GraphChecker gc {graph};
            gc.trace(ENTER);
            // debug(gc.covered);
            debug(1.0 * gc.get_covered_area() / total_cells);
        }

        vector<string> result;

        double total_time = get_time() - start;
        debug(total_time);
        return result;
    }

    static string format_result(int x, int y, char label) {
        ostringstream out;
        out << y << " " << x << " " << label;
        return out.str();
    }
};
