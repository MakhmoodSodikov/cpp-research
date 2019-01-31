#include <iostream>
#include <string>
#include <vector>
#include <algorithm>


class ZFunction {
public:
    ZFunction(std::string& Template, std::string& Text);
    void getTemplatePositions() const;
private:
    const std::string g_sSEPARATOR = "#";
    std::string m_sTemplate;
    std::string m_sText;
};


ZFunction::ZFunction(std::string& Template, std::string& Text) {
    this->m_sTemplate = Template;
    this->m_sText = Text;
}


void ZFunction::getTemplatePositions() const {
//  std::vector<unsigned> positions;
    std::string frame = m_sTemplate + g_sSEPARATOR + m_sText;
    std::vector<unsigned> z_func(frame.length(), 0);
    unsigned left_iter = 0, right_iter = 0;

//  setting z[0] to length of frame
    z_func[0] = frame.length();

    for (unsigned iter = 1; iter < frame.length(); ++iter) {

        if (iter <= right_iter) {
//          we can set z[i] to min(r-l-1, z[i-l]) for initial accuracy of z[i]
//          because we know, that frame[0...r-l] equals to frame[i...z[i]+i-l]
            z_func[iter] = std::min(right_iter - iter + 1, z_func[iter - left_iter]);

//!         UPD: we do not access the memory stack more than O(p), because naive alg works not more than O(p)
//          naive algorithm
            while (iter + z_func[iter] < frame.length() && frame[z_func[iter]] == frame[iter + z_func[iter]]) {
                z_func[iter]++;
            }
        } else {
//!         UPD: we do not access the memory stack more than O(p)
//          naive algorithm if we don't know anything about symbols after iter position in frame
            while (iter + z_func[iter] < frame.length() && frame[z_func[iter]] == frame[iter + z_func[iter]]) {
                z_func[iter]++;
            }
        }

//      reinitializing left and right borders of maximal z-block
        if (iter + z_func[iter] - 1 > right_iter) {
            right_iter = iter + z_func[iter] - 1;
            left_iter = iter;
        }

//!     we don't save all entering positions, so we don't spend more than O(1) memory ¯\_(ツ)_/¯
//      pushing back position into result array if it's value equals to template string length
        if (z_func[iter] == m_sTemplate.length()) {
            std::cout << (iter - m_sTemplate.length() - 1);
        }
    }

//! return positions;
}


//  overload standard output stream
std::ostream& operator << (std::ostream& os, const std::vector<unsigned>& vec) {
    for (auto elem : vec) {
        std::cout << elem << ' ';
    }
    return os;
}


int main() {
    std::string tmp, str;

//  reading data
    std::cin >> tmp;
    std::cin >> str;

//  declaring instance of z-function class
    ZFunction z_function(tmp, str);

//  printing all positions in getTemplatePositions() array
    z_function.getTemplatePositions();

    return 0;
}
