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


const int RED = 31;
const int GREEN = 32;
const int YELLOW = 33;
const int BLUE = 34;
const int MAGENTA = 35;
const int CYAN = 36;
const int WHITE = 37;

void draw_maze(std::function<int (PackedCoord)> color_fn) {
    return;
    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            if (maze[y][x] == 'E')
                cerr << "\033[7";  // inverse
            else
                cerr << "\033[0";
            int color = color_fn(pack(x, y));
            if (color != -1) {
                assert(color >= 30);
                assert(color <= 37);
                cerr << ";" << color;
            }
            cerr << "m";

            cerr << maze[y][x];
            cerr << "\033[0m";  // default
            cerr << " ";
        }
        cerr << endl;
    }
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

    int code() const {
        return pack(x, y) * 4 + dir;
    }
};

std::ostream& operator<<(std::ostream &out, Spin s) {
    out << "Spin(" << s.x << ", " << s.y << ", " << s.dir << ")";
    return out;
}


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


typedef vector<int> Counter;  // indices are PackedCoords


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

    int new_in_counter(const Counter &counter) const {
        int result = 0;
        for (PackedCoord p : cells) {
            if (counter[p] == 0)
                result++;
        }
        return result;
    }

    int add_to_counter(Counter &counter) const {
        ensure_sorted();
        int result = 0;
        for (PackedCoord p : cells) {
            if (counter[p] == 0)
                result++;
            counter[p]++;
        }
        return result;
    }

    void subtract_from_counter(Counter &counter) const {
        for (PackedCoord p : cells) {
            counter[p]--;
        }
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


struct Edge {
private:
    mutable Vertex _from;
    mutable Vertex _to;
    map<PackedCoord, char> _edits;
    Spin _start_spin;
    Spin _end_spin;
    CellSet _path_cells;

public:
    Edge() = default;
    Edge(
        Vertex from,
        Vertex to,
        const map<PackedCoord, char> &edits,
        Spin start_spin,
        Spin end_spin,
        const vector<PackedCoord> path)
        : _from(from), _to(to), _edits(edits),
          _start_spin(start_spin), _end_spin(end_spin)
    {
        for (PackedCoord p : path)
            _path_cells.add_cell(p);
    }

    bool contradicts(const Edge &other) const {
        for (auto kv : _edits)
            if (other._path_cells.contains(kv.first))
                return true;
        for (auto kv : other._edits)
            if (_path_cells.contains(kv.first))
                return true;
        return false;
    }

    int from() const { return _from; }
    int to() const { return _to; }
    int path_length() const { return _path_cells.size(); }
    const CellSet& path_cells() const { return _path_cells; }
    Spin start_spin() const { return _start_spin; }
    Spin end_spin() const { return _end_spin; }
    const map<PackedCoord, char>& edits() const { return _edits; };

    // it's constness makes is a hack
    void flip() const {
        swap(_from, _to);
        if (_from == EXIT)
            _from = ENTER;
        if (_to == ENTER)
            _to = EXIT;
    }
};

ostream& operator<<(ostream &out, const Edge &e) {
    out << "Edge(from=" << e.from()
        << ", to=" << e.to()
        << ", edits=" << e.edits() << ")";
    return out;
}


class EdgeFinder {
public:
    vector<Edge> edges;
private:
    vector<PackedCoord> current_path;
    map<PackedCoord, char> edits;
    Spin start_spin;
    Vertex from;

    void add_edge(Spin end_spin, Vertex to) {
        edges.emplace_back(from, to, edits, start_spin, end_spin, current_path);
    }

    void trace(Spin s, int num_edits) {
        char c = maze[s.y][s.x];
        if (c == '.' || c == 'E') {
            if (num_edits == 0) {
                Vertex to = c == '.' ? EXIT : pack(s.x, s.y);
                if (from != to)
                    add_edge(s, to);
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
                if (d == new_dir)
                    continue;
                if (current_path.size() > 1 && (d ^ s.dir) == 2)
                    continue;
                edits[pos] = "SLUR"[(4 + d - s.dir) % 4];
                new_spin = s;
                new_spin.dir = d;
                new_spin.advance();
                trace(new_spin, num_edits - 1);
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


class GraphChecker {
private:
    map<Vertex, vector<const Edge*>> graph;
    vector<Vertex> current_path;
    vector<const Edge*> current_edges;
public:
    CellSet covered;
    set<PackedCoord> forward_reachable;

    GraphChecker(const vector<const Edge*> edges) {
        for (auto e : edges) {
            graph[e->from()].push_back(e);
        }
    }

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

        forward_reachable.insert(v);

        auto kv = graph.find(v);
        if (kv == graph.end()) {
            return;
        }

        for (auto e : kv->second) {
            Vertex w = e->to();
            if (find(current_path.begin() + 1, current_path.end(), w) !=
                current_path.end())
                continue;
            bool overlaps = false;
            for (const Edge* pe : current_edges) {
                if (e->path_cells().overlaps(pe->path_cells())) {
                    overlaps = true;
                    break;
                }
            }
            if (overlaps)
                continue;

            current_path.push_back(w);
            current_edges.push_back(e);

            rec();

            assert(current_path.back() == w);
            current_path.pop_back();
            assert(current_edges.back() == e);
            current_edges.pop_back();
        }
    }
};


vector<Edge> candidate_edges;


typedef vector<const Edge*> Graph;


set<Vertex> reachable(const Graph &graph, bool forward) {
    set<Vertex> result;
    if (forward)
        result.insert(ENTER);
    else
        result.insert(EXIT);

    while (true) {
        bool changed = false;
        for (auto e : graph) {
            int from = forward ? e->from() : e->to();
            int to = forward ? e->to() : e->from();
            if (result.count(from) > 0 && result.count(to) == 0) {
                result.insert(to);
                changed = true;
            }
        }
        if (!changed)
            break;
    }
    return result;
}

set<Vertex> compute_forward_reachable(const Graph &graph) {
    GraphChecker gc {graph};
    gc.trace(ENTER);
    return gc.forward_reachable;
}

set<Vertex> compute_backward_reachable(const Graph &graph) {
    for (auto e : graph)
        e->flip();
    GraphChecker gc {graph};
    gc.trace(ENTER);
    for (auto e : graph)
        e->flip();
    return gc.forward_reachable;
}


Graph greedy_fill(const Graph &other_edges, PackedCoord pt, int budget, bool prune_unreachable) {
    assert(maze[unpack_y(pt)][unpack_x(pt)] == 'E');
    if (budget > 19)
        budget = 19;
    vector<Graph> incoming_candidates(20);
    vector<Graph> outgoing_candidates(20);

    for (const Edge &e : candidate_edges) {
        if (e.from() != pt && e.to() != pt)
            continue;

        bool contradicts = false;
        for (auto e2 : other_edges)
            if (e.contradicts(*e2)) {
                contradicts = true;
                break;
            }
        if (contradicts)
           continue;

        if (e.from() == pt)
            outgoing_candidates[e.edits().size()].push_back(&e);
        if (e.to() == pt)
            incoming_candidates[e.edits().size()].push_back(&e);
    }

    /*for (auto candidates : {&incoming_candidates, &outgoing_candidates}) {
        for (auto &ic : *candidates) {
            sort(ic.begin(), ic.end(), [](const Edge *e1, const Edge *e2){
                return e1->path_cells().size() > e2->path_cells().size();
            });
            if (ic.size() > 20)
                ic.resize(20);
        }
    }*/

    set<Vertex> forward_reachable;
    if (prune_unreachable)
        forward_reachable = compute_forward_reachable(other_edges);
    set<Vertex> backward_reachable;
    if (prune_unreachable)
        backward_reachable = compute_backward_reachable(other_edges);
    // debug(forward_reachable);
    // debug(backward_reachable);

    // debug(outgoing_candidates.size());
    // debug(incoming_candidates.size());

    Graph best_solution;
    double best_solution_score = -100;

    Counter counter(pack(W - 1, H - 1));
    for (auto oe : other_edges)
        oe->path_cells().add_to_counter(counter);

    Graph incoming;
    Graph outgoing;
    int covered = 0;

    auto contradicts_existing= [&](const Edge *e) {
        for (auto e2 : incoming) {
            if (e == e2 || e->contradicts(*e2))
                return true;
        }
        for (auto e2 : outgoing) {
            if (e == e2 || e->contradicts(*e2))
                return true;
        }
        return false;
    };

    std::function<void(int)> rec;
    rec = [&](int budget) {
        assert(budget >= 0);

        if (!outgoing.empty() && !incoming.empty()) {
            // TODO: more accurate counting
            set<Vertex> endpoints;

            double score = covered;
            for (auto e : incoming) {
                //score += e->path_cells().size();
                endpoints.insert(e->from());

                if (forward_reachable.count(e->from()) > 0)
                    score += 10;
            }
            for (auto e : outgoing) {
                //score += e->path_cells().size();
                endpoints.insert(e->to());

                if (backward_reachable.count(e->from()) > 0)
                    score += 10;
            }
            score += 0.001 * budget;
            score += 0.1 * endpoints.size();

            if (score > best_solution_score) {
                best_solution_score = score;
                best_solution = incoming;
                copy(outgoing.begin(), outgoing.end(),
                     back_inserter(best_solution));
            }
        }

        if (outgoing.size() < incoming.size() + 0.01 * (rand() % 2)) {
            // adding outgoing

            for (int i = 0; i <= budget; i++) {
                const Edge *best_edge = nullptr;
                double best_score = -1000;
                for (auto new_edge : outgoing_candidates[i]) {
                    // if (prune_unreachable &&
                    //     backward_reachable.count(new_edge->to()) == 0)
                    //     continue;

                    if (contradicts_existing(new_edge))
                        continue;

                    bool reachable = incoming.size() == 0;
                    for (auto ie : incoming) {
                        if (ie->from() != new_edge->to() &&
                            !ie->path_cells().overlaps(new_edge->path_cells())) {
                            reachable = true;
                            break;
                        }
                    }
                    if (!reachable)
                        continue;

                    double score = new_edge->path_cells().new_in_counter(counter);
                    if (score > best_score) {
                        best_score = score;
                        best_edge = new_edge;
                    }
                }

                if (best_edge != nullptr) {
                    outgoing.push_back(best_edge);
                    int delta_covered = best_edge->path_cells().add_to_counter(counter);
                    covered += delta_covered;

                    rec(budget - best_edge->edits().size());

                    assert(outgoing.back() == best_edge);
                    outgoing.pop_back();
                    best_edge->path_cells().subtract_from_counter(counter);
                    covered -= delta_covered;
                }
            }
        } else {
            // adding incoming

            for (int i = 0; i <= budget; i++) {
                const Edge *best_edge = nullptr;
                double best_score = -1000;
                for (auto new_edge : incoming_candidates[i]) {
                    // if (prune_unreachable &&
                    //     forward_reachable.count(new_edge->from()) == 0)
                    //     continue;

                    if (contradicts_existing(new_edge))
                        continue;

                    bool reachable = outgoing.size() == 0;
                    for (auto oe : outgoing) {
                        if (oe->to() != new_edge->from() &&
                            !oe->path_cells().overlaps(new_edge->path_cells())) {
                            reachable = true;
                            break;
                        }
                    }
                    if (!reachable)
                        continue;

                    double score = new_edge->path_cells().new_in_counter(counter);
                    if (score > best_score) {
                        best_score = score;
                        best_edge = new_edge;
                    }
                }

                if (best_edge != nullptr) {
                    incoming.push_back(best_edge);
                    int delta_covered = best_edge->path_cells().add_to_counter(counter);
                    covered += delta_covered;

                    rec(budget - best_edge->edits().size());

                    assert(incoming.back() == best_edge);
                    incoming.pop_back();
                    best_edge->path_cells().subtract_from_counter(counter);
                    covered -= delta_covered;
                }
            }
        }
    };

    rec(budget);
    assert(covered == 0);

/*
    debug(best_solution_score);
    debug(best_solution.size());
    for (auto e : best_solution)
        debug(*e);

    draw_maze([&](PackedCoord p) {
        if (p == pt)
            return CYAN;
        bool out = false;
        bool in = false;

        for (auto e : best_solution) {
            if (e->edits().count(p) > 0)
                return BLUE;
            if (e->to() == pt) {
                if (e->path_cells().contains(p) || e->from() == p)
                    in = true;
            } else {
                assert(e->from() == pt);
                if (e->path_cells().contains(p) || e->to() == p)
                    out = true;
            }
        }

        if (in && out)
            return YELLOW;
        if (in && !out)
            return GREEN;
        if (!in && out)
            return RED;
        return -1;
    });*/

    return best_solution;
}

void remove_vertex(Graph &graph, PackedCoord pt) {
    auto end = remove_if(graph.begin(), graph.end(), [pt](const Edge *e){
        return e->from() == pt || e->to() == pt;
    });
    graph.erase(end, graph.end());
}


class MazeFixing {
public:
    vector<string> improve(vector<string> maze, int F) {
        ::H = maze.size();
        ::W = maze.front().size();
        ::F = F;
        ::maze = maze;
        ::original_maze = maze;

        double start_time = get_time();

        debug3(W, H, F);

        int total_cells = 0;
        for (const string &row : maze)
            total_cells += row.size() - count(row.begin(), row.end(), '.');

        vector<PackedCoord> es;
        for (int y = 0; y < H; y++) {
            for (int x = 0; x < W; x++) {
                if (maze[y][x] == 'E')
                    es.push_back(pack(x, y));
            }
        }
        debug(es.size());

        int qq = 0;

        int total_candidates_size = 0;

        candidate_edges.clear();
        for (int num_edits = 0; num_edits <= 5; num_edits++) {
            if (num_edits == 5 && es.size() > 30)
                break;
            debug(num_edits);

            map<pair<int, int>, int> cnt;
            map<int, int> lengths;

            EdgeFinder ef;

            for (auto p : enumerate_starting_points()) {
                Vertex start = p.first;
                ef.trace(start, p.second, num_edits);

                if (++qq % 50 == 0 && get_time() > start_time + TIME_LIMIT * 0.4) {
                    cerr << "candidates time limit" << endl;
                    break;
                }
            }
            for (Edge &e : ef.edges) {
                candidate_edges.push_back(e);
                total_candidates_size += e.path_cells().size();

                cnt[{e.from(), e.to() == EXIT}]++;
                lengths[e.path_length()]++;
            }
            debug(cnt);
            // debug(cnt.size());
            // debug(lengths);
        }
        debug(candidate_edges.size());
        debug(total_candidates_size);

        vector<const Edge*> actual_edges;
        for (const Edge &e : candidate_edges) {
            if (e.edits().empty())
                actual_edges.push_back(&e);
        }
        GraphChecker gc {actual_edges};
        gc.trace(ENTER);
        // debug(gc.covered);

        double predicted_score = 1.0 * gc.get_covered_area() / total_cells;
        debug(predicted_score);

        /*draw_maze([&](PackedCoord p) {
            if (gc.covered.contains(p))
                return RED;
            else if (gc.forward_reachable.count(p) > 0)
                return GREEN;
            return -1;
        });*/

        if (W * H > 1500) {
            vector<string> result;
            /// Greedy solution

            sort(candidate_edges.begin(), candidate_edges.end(),
                [](const Edge &e1, const Edge &e2){
                    return e1.path_cells().size() > e2.path_cells().size();
                });
            int remaining_edits = F;
            Graph solution;
            for (const Edge &e : candidate_edges) {
                if (e.from() > e.to() &&
                    e.from() != ENTER &&
                    e.to() != EXIT)
                    continue;

                bool contradicts = false;
                for (const Edge *pe : solution) {
                    if (pe->contradicts(e)) {
                        contradicts = true;
                        break;
                    }
                }
                if (contradicts)
                    continue;

                if (e.edits().size() <= remaining_edits) {
                    solution.push_back(&e);
                    remaining_edits -= e.edits().size();

                    if (get_time() > start_time + TIME_LIMIT) {
                        cerr << "TIME LIMIT" << endl;
                        break;
                    }
                }
            }

            for (const Edge *e : solution) {
                for (auto kv : e->edits()) {
                    PackedCoord p = kv.first;
                    result.push_back(
                        format_result(unpack_x(p), unpack_y(p), kv.second));
                }
            }
            draw_maze([&](PackedCoord p) {
                for (const Edge *e : solution)
                    if (e->path_cells().contains(p))
                        return RED;
                return -1;
            });

            debug2(result.size(), F);

            double total_time = get_time() - start_time;
            debug(total_time);

            return result;
        }

        vector<string> result;

        Graph solution;
        for (int i = 0; i < 2; i++) {
            debug(i);
            for (PackedCoord pt : es) {
                remove_vertex(solution, pt);

                int budget = F;
                for (auto e : solution)
                    budget -= e->edits().size();

                auto new_edges = greedy_fill(
                    solution, pt,
                    min<int>(2 * F / es.size(), budget),
                    /*prune_unreachable*/ i > 0);
                copy(new_edges.begin(), new_edges.end(),
                     back_inserter(solution));

                if (++qq % 30 == 0 && get_time() > start_time + TIME_LIMIT * 0.8) {
                    cerr << "greedy time limit" << endl;
                    break;
                }
            }

            set<Vertex> forward_reachable = compute_forward_reachable(solution);
            set<Vertex> backward_reachable = compute_backward_reachable(solution);

            draw_maze([&](PackedCoord p) {
                for (auto e : solution) {
                    if (e->path_cells().contains(p))
                        return BLUE;
                }

                bool out = forward_reachable.count(p) > 0;
                bool in = backward_reachable.count(p) > 0;
                if (in && out)
                    return YELLOW;
                if (in && !out)
                    return GREEN;
                if (!in && out)
                    return RED;
                return -1;
            });

        }

        debug(solution.size());

        Counter counter(pack(W - 1, H - 1));
        for (auto oe : solution)
            oe->path_cells().add_to_counter(counter);

        auto contradicts_existing = [&](const Edge *e) {
            for (auto e2 : solution) {
                if (e == e2 || e->contradicts(*e2))
                    return true;
            }
            return false;
        };

        set<Vertex> forward_reachable = compute_forward_reachable(solution);
        set<Vertex> backward_reachable = compute_backward_reachable(solution);
        while (true) {
            if (++qq % 2 == 0 && get_time() > start_time + TIME_LIMIT) {
                cerr << "improvement time limit" << endl;
                break;
            }

            int budget = F;
            for (auto e : solution)
                budget -= e->edits().size();

            double best_score = 0;
            const Edge *best_edge = nullptr;
            for (const auto &e : candidate_edges) {
                if (e.edits().size() > budget)
                    continue;
                if (e.path_cells().size() / (e.edits().size() + 2.0) <= best_score)
                    continue;
                if (contradicts_existing(&e))
                    continue;
                if (forward_reachable.count(e.from()) == 0)
                    continue;
                if (backward_reachable.count(e.to()) == 0)
                    continue;

                double score = e.path_cells().new_in_counter(counter) / (e.edits().size() + 2.0);
                if (score > best_score) {
                    best_score = score;
                    best_edge = &e;
                }
            }

            if (best_edge == nullptr)
                break;

            debug2(best_score, *best_edge);
            solution.push_back(best_edge);
            best_edge->path_cells().add_to_counter(counter);
        }


        for (const Edge *e : solution) {
            for (auto kv : e->edits()) {
                PackedCoord p = kv.first;
                result.push_back(
                    format_result(unpack_x(p), unpack_y(p), kv.second));
            }
        }
        draw_maze([&](PackedCoord p) {
            for (const Edge *e : solution)
                if (e->path_cells().contains(p))
                    return GREEN;
            return -1;
        });


        debug2(result.size(), F);

        double total_time = get_time() - start_time;
        debug(total_time);

        return result;
    }

    static string format_result(int x, int y, char label) {
        ostringstream out;
        out << y << " " << x << " " << label;
        return out.str();
    }
};
