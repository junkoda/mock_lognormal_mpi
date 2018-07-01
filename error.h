#ifndef ERROR_H
#define ERROR_H 1

class Error{};

class IOError : public Error {};
class FileNotFoundError : public IOError{};
class TypeError : public Error{};
class ValueError : public Error{};
class RuntimeError : public Error{};

#endif
