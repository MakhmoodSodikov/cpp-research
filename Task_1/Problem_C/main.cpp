#include <iostream>
#include <vector>


char const Q_SIGN = '?';
int const NIL = -1;
int const ALPHABET_SIZE = 26;
char const ALPHA = 'a';


char get_char_code (char x) {
    return x - ALPHA;
}


std::ostream& operator << (std::ostream& os, const std::vector<int>& vec) {
    for (auto elem : vec) {
        std::cout << elem << ' ';
    }
    return os;
}


struct Vertex {
    Vertex(int parent, char character)
    : m_nParentChar(character),
    m_nParent(parent),
    m_nSuffixLink(NIL),
    m_bIsPattern(false) {}

    std::vector<int> m_vEdges = std::vector<int>(ALPHABET_SIZE, NIL);
    std::vector<int> m_vFastMove = std::vector<int>(ALPHABET_SIZE, NIL);
    std::vector<int> m_vPatternPos;
    int m_nSuffixLink;
    int m_nParent;
    char m_nParentChar;
    bool m_bIsPattern;
};


class AhoCorasickAutomaton {
public:
    explicit AhoCorasickAutomaton(const std::string &mask);
    friend class FindMatches;
    void addPattern(const std::pair<int, int> &submask_pos, int patt_index);
    void findSubmaskPositions(const std::string &mask);
    int getSuffLink(int index);
    int go(int index, char character);
    void addPatterns();
    void checkFirstSymbols(std::pair<int, int> &pos, const std::string &mask);
    void checkCentralSymbols(std::pair<int, int> &pos, const std::string &mask);
    void checkLastSymbols(std::pair<int, int> &pos, const std::string &mask);
private:
    std::string m_sMask;
    std::vector<Vertex> m_Trie;
    std::vector<std::pair<int, int>> m_vSubMasPos;
    void checkNIL(int &vert, char &character);
};


class FindMatches {
private:
    std::string m_sText;
    AhoCorasickAutomaton m_Aho;
    std::vector<int> result;
    void findEntrings();
    void increaseCounters(int i, int u, std::vector<int>& entries);
    void precountResults(std::vector<int>& entries);
public:
    explicit FindMatches(std::string& mask, std::string& text) : m_Aho(mask), m_sText(text) {};
    void solve();
};


AhoCorasickAutomaton::AhoCorasickAutomaton(const std::string &mask)
    : m_Trie(1, Vertex(0, NIL)),
    m_sMask(mask) {

    //initializing root-node
    m_Trie[0].m_nSuffixLink = 0;

    findSubmaskPositions(mask);
    addPatterns();
}


void AhoCorasickAutomaton::addPattern(const std::pair<int, int> &submask_pos, int patt_index) {
    int curr_vert = 0;

    for (int i = submask_pos.first; i <= submask_pos.second; i++) {
        char character = get_char_code(m_sMask[i]);
        checkNIL(curr_vert, character);

        curr_vert = m_Trie[curr_vert].m_vEdges[character];
    }

    m_Trie[curr_vert].m_bIsPattern = true;

    m_Trie[curr_vert].m_vPatternPos.push_back(patt_index);

}


int AhoCorasickAutomaton::getSuffLink(int index) {
    if (m_Trie[index].m_nSuffixLink == NIL) {

        if (m_Trie[index].m_nParent == 0) {

            m_Trie[index].m_nSuffixLink = 0;
        } else {
            auto suff_link = getSuffLink(m_Trie[index].m_nParent);

            m_Trie[index].m_nSuffixLink = go(suff_link, m_Trie[index].m_nParentChar);

        }
    }

    return m_Trie[index].m_nSuffixLink;

}


void AhoCorasickAutomaton::findSubmaskPositions(const std::string &mask) {
    std::pair<int, int> current_submask_pos;

    checkFirstSymbols(current_submask_pos, mask);

    checkCentralSymbols(current_submask_pos, mask);

    checkLastSymbols(current_submask_pos, mask);

}


int AhoCorasickAutomaton::go(int index, char character) {
    if (m_Trie[index].m_vFastMove[character] == NIL) {

        if (m_Trie[index].m_vEdges[character] != NIL) {

            m_Trie[index].m_vFastMove[character] = m_Trie[index].m_vEdges[character];

        } else if (index == 0) {

            m_Trie[index].m_vFastMove[character] = 0;

        } else {

            m_Trie[index].m_vFastMove[character] = go(getSuffLink(index), character);

        }
    }

    return m_Trie[index].m_vFastMove[character];
}


void AhoCorasickAutomaton::addPatterns() {
    for (int i = 0; i < m_vSubMasPos.size(); i++) {

        addPattern(m_vSubMasPos[i], i);

    }
}


void AhoCorasickAutomaton::checkFirstSymbols(std::pair<int, int>& pos, const std::string &mask) {
    if (isalpha(mask[0])) {

        pos.first = 0;

    }

    if (m_sMask[1] == Q_SIGN && isalpha(m_sMask[0])) {

        pos.second = 0;
        m_vSubMasPos.push_back(pos);

    }

}


void AhoCorasickAutomaton::checkCentralSymbols(std::pair<int, int>& pos, const std::string &mask) {
    for (int i = 1; i < m_sMask.length() - 1; i++) {
        if (mask[i - 1] == Q_SIGN && isalpha(mask[i])) {

            pos.first = i;

        }

        if (mask[i + 1] == Q_SIGN && isalpha(mask[i])) {

            pos.second = i;
            m_vSubMasPos.push_back(pos);

        }
    }
}


void AhoCorasickAutomaton::checkLastSymbols(std::pair<int, int>& pos, const std::string &mask) {

    if (m_sMask[mask.length() - 2] == Q_SIGN && isalpha(mask[mask.length() - 1])) {

        pos.first = (int)mask.length() - 1;

    }

    if (isalpha(m_sMask[m_sMask.length() - 1])) {

        pos.second = (int)(mask.length() - 1);
        m_vSubMasPos.push_back(pos);

    }

}


void AhoCorasickAutomaton::checkNIL(int &vert, char &character) {
    if (m_Trie[vert].m_vEdges[character] == NIL) {

        m_Trie.emplace_back(Vertex(vert, character));

        m_Trie[vert].m_vEdges[character] = (int)(m_Trie.size() - 1);

    }
}


void FindMatches::solve() {
    findEntrings();
    std::cout << result;
}


void FindMatches::findEntrings() {
    std::vector<int> entries(m_sText.length());

    int v = 0;

    for (int i = 0; i < m_sText.length(); i++) {
        v = m_Aho.go(v, get_char_code(m_sText[i]));
        int u = v;

        do {
            if (m_Aho.m_Trie[u].m_bIsPattern) {
                increaseCounters(i, u, entries);
            }

            u = m_Aho.getSuffLink(u);
        } while (u != 0);

    }

    precountResults(entries);
}


void FindMatches::increaseCounters(int i, int u, std::vector<int>& entries) {
    for (int index : m_Aho.m_Trie[u].m_vPatternPos) {
        int b = m_Aho.m_vSubMasPos[index].second;
        int a = m_Aho.m_vSubMasPos[index].first;
        int len = m_Aho.m_sMask.length();

        int startIndex = i - b + a;

        if ((startIndex - a >= 0) && (startIndex - a + len - 1 < m_sText.length())) {
            entries[startIndex - a]++;
        }
    }
}


void FindMatches::precountResults(std::vector<int>& entries) {
    for (int i = 0; i < entries.size(); i++) {

        if (entries[i] == m_Aho.m_vSubMasPos.size()) {
            result.push_back(i);
        }

    }
}


int main() {
    std::string pattern, text;
    std::cin >> pattern >> text;
    FindMatches fm(pattern, text);
    fm.solve();
    return 0;
}

/*
 * tests:
 *
a?aa?aab
aaaaaaaabaab

?a?
aaa

a?a?a?a?a
abaaaaaaaab

ab??
dcabd

??ab
dabd

ab??aba
ababacaba
 */