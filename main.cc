#include <iostream>
#include <vector>
#include <string>
#include <cassert>

using namespace std;


#include "solution.h"


int main(int argc, char **argv) {
    (void)argc; (void)argv;  // suppress unused parameter warning

    int H;
    cin >> H;

    vector<string> maze(H);
    for (auto &row : maze)
        cin >> row;

    int F;
    cin >> F;

    auto result = MazeFixing().improve(maze, F);

    cout << result.size() << endl;
    for (auto s : result)
        cout << s.c_str() << endl;

    cout.flush();
    return 0;
}
