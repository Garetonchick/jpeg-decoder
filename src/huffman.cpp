#include <huffman.h>
#include <cstdint>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>

class HuffmanTree::Impl {
public:
    void Build(const std::vector<uint8_t> &code_lengths, const std::vector<uint8_t> &values) {
        if (code_lengths.size() > 16) {
            throw std::invalid_argument("HuffmanTree::Build: code_lenghts is too large");
        }

        if (std::accumulate(code_lengths.begin(), code_lengths.end(), size_t{}) != values.size()) {
            throw std::invalid_argument(
                "HuffmanTree::Build: values size doesn't match code_lengths sum");
        }

        Clear();

        uint16_t code = 0;
        uint16_t val_idx = 0;

        for (uint8_t len = 1; len <= code_lengths.size(); ++len) {
            for (uint8_t i = 0; i < code_lengths[len - 1]; ++i) {
                Insert(values[val_idx], code, len);
                ++code;
                ++val_idx;
            }

            code <<= 1;
        }
    }

    bool Move(bool bit, int &value) {
        if (tree_.at(current_node_).nxt[bit] == kNullPtr) {
            throw std::invalid_argument("HuffmanTree::Move: Invalid move");
        }

        current_node_ = tree_[current_node_].nxt[bit];

        if (tree_[current_node_].has_value) {
            value = tree_[current_node_].val;
            current_node_ = 0;
            return true;
        }

        return false;
    }

private:
    void Insert(uint8_t val, uint16_t bit_path, uint8_t len) {
        uint32_t cur = 0;

        while (len--) {
            uint8_t bit = ((bit_path >> len) & 1);

            if (tree_[cur].nxt[bit] == kNullPtr) {
                tree_[cur].nxt[bit] = NewNode();
            }

            cur = tree_[cur].nxt[bit];

            if (tree_[cur].has_value) {
                throw std::invalid_argument(
                    "HuffArchiver::Insert: Met terminal node while inserting");
            }
        }

        tree_[cur].val = val;
        tree_[cur].has_value = true;
    }

    uint32_t NewNode() {
        tree_.push_back({});
        return tree_.size() - 1;
    }

    void Clear() {
        tree_.resize(1);
        tree_[0] = TreeNode{};
        current_node_ = 0;
    }

private:
    static const uint32_t kNullPtr = std::numeric_limits<uint32_t>::max();

    struct TreeNode {
        uint32_t nxt[2] = {kNullPtr, kNullPtr};
        uint8_t val = 0;
        bool has_value = false;
    };

private:
    std::vector<TreeNode> tree_ = {TreeNode{}};
    uint32_t current_node_ = 0;
};

HuffmanTree::HuffmanTree() {
    impl_ = std::make_unique<HuffmanTree::Impl>();
}

void HuffmanTree::Build(const std::vector<uint8_t> &code_lengths,
                        const std::vector<uint8_t> &values) {
    impl_->Build(code_lengths, values);
}

bool HuffmanTree::Move(bool bit, int &value) {
    return impl_->Move(bit, value);
}

HuffmanTree::HuffmanTree(HuffmanTree &&) = default;

HuffmanTree &HuffmanTree::operator=(HuffmanTree &&) = default;

HuffmanTree::~HuffmanTree() = default;
