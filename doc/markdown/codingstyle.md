# C++ Coding Style Guide {#codingstyle}

[TOC]

## 1. Introduction {#codingstyle-introduction}
This document outlines the coding standards and best practices for writing Geneasic's C++ code. Consistent coding style helps improve code readability, maintainability, and reduces the likelihood of bugs. Please follow the guides bellow.

## 2. File Organization {#codingstyle-file-organization}

### 2.1. File Naming
* Use lowercase letters.
* Separate words with underscores (_).
* Header files: .hpp.
* Source files: .cpp.

**Example:**
```cpp
my_class.hpp
my_class.cpp
```

## 3. Naming Conventions {#codingstyle-naming-conventions}

### 3.1. Variables
* Use meaningful names.
* Separate words with underscores (_).

**Example:**

```cpp
int user_age;
```

### 3.2 Functions
* Separate words with underscores (_).
* Use `auto` if the name of return type is too long.
    * Optional:
        * Could use `auto` with `->` if the name of return type is too long.

**Example:**
```cpp
int calculate_sum(const int &a, const &b) {
  return a + b;
};

// Using auto
auto do_something(const int &a, const int &b, const int &c) {
  std::unordered_map<int, std::vector<std::vector<double>>> res;
  for (int i=0; i!=c; ++i) {
    double num = i;
      res[i] = std::vector<std::vector<double>>(a, std::vector<double>(b, 0));
    for (auto &values: res[i])
      for (auto &value: values)
        value = num *= i;
  }
  return res;
}

// Accetable
auto do_something2(const int &a, const int &b, const int &c) 
  -> std::unordered_map<int, std::vector<std::vector<double>>> {
  std::unordered_map<int, std::vector<std::vector<double>>> res;
  for (int i=0; i!=c; ++i) {
    double num = i;
      res[i] = std::vector<std::vector<double>>(a, std::vector<double>(b, 0));
    for (auto& values: res[i])
      for (auto& value: values)
        value = num *= i;
  }
  return res;
}
```

### 3.3. Classes and Structs

* Use PascalCase.

**Example:**
```cpp
class MyClass;
```
### 3.4 Constants 
* Use ALL_CAPS with underscores for constants.
* Separate words with underscores (_).

**Example:**
```cpp
static constexpr int MAX_USERS = 100;
```

## 4. Formatting {#codingstyle-formatting}

### 4.1. Indentation
* Use 2 spaces per indentation level. Do not use tabs.

### 4.2. Braces

* Use K&R style
    * In this style, a function has its opening and closing braces on their own lines and with the same indentation as the declaration, and with statements indented an additional level than the declaration. A multi-statement block inside a function, however, has its opening brace on the same line as its control clause while the closing brace remains on its own line unless followed by a keyword such as `else`.

**Example:**
```cpp
int main(int argc, char *argv[]) {
  while (x == y) {
    do_something();
    if (some_error)
      fix_issue(); // single-statement block without braces
    else {
      continue_as_usual();
      do_something_else();
    }
  }
  auto done = final_thing();
    
  // New line for else (or else if) is acceptable
  if (done) {
    celebrate();
    std::cout << "Sucess\n";
  }
  else {
    crying();
    std::cerr << "Error\n";
  }
  return 0;
}
```

### 4.3. Line Length

* Limit lines to 100 characters.


### 4.4. Whitespace

* Use spaces around operators.
* Do not add spaces after opening or before closing parentheses in function calls.

**Example:**
```cpp
int sum = a + b;
void myFunction(int x, int y);
```

## 5. Comments {#codingstyle-comments}

### 5.1. Using doxygen's comment style 

* [Doxygen](https://github.com/doxygen/doxygen)

**Example:**

```cpp
/**
 *  A test class. A more elaborate class description.
 */
 
class Javadoc_Test {

public:
    
  /**
   * A constructor.
   * A more elaborate description of the constructor.
   */
  Javadoc_Test();
 
  /**
   * A destructor.
   * A more elaborate description of the destructor.
   */
 ~Javadoc_Test();
    
  /**
   * a normal member taking two arguments and returning an integer value.
   * @param a an integer argument.
   * @param s a constant character pointer.
   * @see Javadoc_Test()
   * @see ~Javadoc_Test()
   * @see testMeToo()
   * @see publicVar()
   * @return The test results
   */
  int testMe(int a,const char *s);
  
  int publicVar; //!< a public variable.
};
```

## 6. Classes and Structs {#codingstyle-classes-and-structs}

### 6.1. Class Definition

* **Private variable should end with underscore.**

**Example:**
```cpp
class MyClass {
public:
  MyClass() = default;
  MyClass(int a): privateVariable1_(a) {}; 
 ~MyClass();
    
  void public_method();

protected:
  void protected_method();

private:
  int privateVariable1_;
  void private_ethod();
};
```

## 7. Templates {#codingstyle-templates}

### 7.1 Template Declaration

* Use descriptive names for template parameters.
* Prefer single uppercase letters for simple, generic types (T, U, V), and more descriptive names for more complex types.

**Example:**
```cpp
template <class T>
class MyClass;

template <class KeyType, class ValueType>
class MyMap;
```

### 7.2 Template Function Definitions

* Define template functions in the same header file, typically at the end or within the class definition.

```cpp
template <class T>
class MyClass {
public:
  void doSomething() {
    // Function definition
  }
};
```

## 8. Miscellaneous {#codingstyle-miscellaneous}

### 8.1 Prefer using standard containers or standard algorithms

```cpp
std::vector<int> numbers{3,4,1,2,3};
std::ranges::sort(numbers); /** 1,2,3,4,5 */
```

### 8.2 Others

* Use "Almost Always Auto" (AAA) for type inference.
* Prefer using the Standard Template Library (STL) over complex nested raw loops.
* Declare member functions of stateless classes/structs as static.
* Utilize const whenever possible.
* Use move semantics instead of copying.
```