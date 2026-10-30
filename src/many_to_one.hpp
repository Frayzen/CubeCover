#include <unordered_map>
#include <unordered_set>

template <typename Key, typename Value, typename KeyHash = std::hash<Key>,
          typename ValueHash = std::hash<Value>>
class ManyToOne {
  // Key -> Value
  std::unordered_map<Key, Value, KeyHash> forward;

  // Value -> {Key1, Key2, ...}
  std::unordered_map<Value, std::unordered_set<Key, KeyHash>, ValueHash>
      reverse;

public:
  void set(const Key &key, const Value &value) {
    remove(key);

    forward[key] = value;
    reverse[value].insert(key);
  }

  void remove(const Key &key) {
    auto it = forward.find(key);
    if (it == forward.end())
      return;

    auto rit = reverse.find(it->second);

    rit->second.erase(key);

    if (rit->second.empty())
      reverse.erase(rit);

    forward.erase(it);
  }

  bool contains(const Key &key) const { return forward.contains(key); }

  const Value &operator[](const Key &key) const { return forward.at(key); }

  const std::unordered_set<Key, KeyHash> &getAll(const Value &value) const {
    return reverse.at(value);
  }

  bool containsValue(const Value &value) const {
    return reverse.contains(value);
  }

  size_t size() const { return forward.size(); }

  void clear() {
    forward.clear();
    reverse.clear();
  }
};
