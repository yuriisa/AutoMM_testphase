#include <exception>
#include <iostream>


namespace exc
{
  class exception : public std::exception
    {
    public:
    exception(int code,
	      const char* file,
	      const char* line,
	      const char* reason);
    virtual void PrintWhat(std::ostream& where) const;
    virtual void PrintWhat() const;
    virtual const char* what() const throw();
    void SetFile(const std::string& file);
    void SetLine(const std::string& line);
    void SetReason(const std::string& reason);
    const std::string& GetReason() const;
    virtual ~exception() throw() {};
    
    protected:
    std::string exc_type;
    
    private:
    int m_code;
    std::string m_file;
    std::string m_line;
    std::string m_reason;
    };

  class bad_alloc : public exception 
  {
  public:
    bad_alloc(int code,
	      const char* file,
	      const char* line,
	      const char* reason) throw() : exception(code, file, line, reason) { exc_type="Bad allocation";}
    virtual ~bad_alloc() throw() {};
  };
  
 class logic_error : public exception 
  {
  public:
    logic_error(int code,
	      const char* file,
	      const char* line,
	      const char* reason) throw() : exception(code, file, line, reason) {exc_type="Logical error";}
    virtual ~logic_error() throw() {};
  };
  
  class invalid_argument : public logic_error
  {
  public:
    invalid_argument(int code,
	      const char* file,
	      const char* line,
		     const char* reason) throw() : logic_error(code, file, line, reason) {exc_type="Invalid argument"; }
    virtual ~invalid_argument() throw() {};
  };

  class domain_error : public logic_error
  {
  public:
    domain_error(int code,
	      const char* file,
	      const char* line,
		 const char* reason) throw() : logic_error(code, file, line, reason) {exc_type="Domain error"; }
    virtual ~domain_error() throw() {};
  };

  class runtime_error : public exception
  {
  public:
    runtime_error(int code,
	      const char* file,
	      const char* line,
	      const char* reason) throw() : exception(code, file, line, reason) {exc_type="Runtime error";  }
    virtual ~runtime_error() throw() {};
  };

  class overflow_error : public runtime_error
  {
  public:
    overflow_error(int code,
	      const char* file,
	      const char* line,
	      const char* reason) throw() : runtime_error(code, file, line, reason) {exc_type="Overflow error";  }
    virtual ~overflow_error() throw() {};
  };

  class range_error : public runtime_error
  {
  public:
    range_error(int code,
	      const char* file,
	      const char* line,
	      const char* reason) throw() : runtime_error(code, file, line, reason) {exc_type="Range error"; }
    virtual ~range_error() throw() {};
  };

  class underflow_error : public runtime_error
  {
  public:
    underflow_error(int code,
	      const char* file,
	      const char* line,
	      const char* reason) throw() : runtime_error(code, file, line, reason) {exc_type="Underflow error"; }
    virtual ~underflow_error() throw() {};
  };

}

#ifndef STRINGIZE
#define STRINGIZE(N) _STRINGIZE(N)

#define _STRINGIZE(N) #N

#endif
#define THROW(code, reason) throw exc::exception(code, __FILE__, STRINGIZE(__LINE__), reason)
#define THROW_BAD_ALLOC(code, reason) throw exc::bad_alloc(code, __FILE__, STRINGIZE(__LINE__), reason)
#define THROW_LOGIC_ERROR(code, reason) throw exc::logic_error(code, __FILE__, STRINGIZE(__LINE__), reason)
#define THROW_INVALID_ARGUMENT(code, reason) throw exc::invalid_argument(code, __FILE__, STRINGIZE(__LINE__), reason)
#define THROW_DOMAIN_ERROR(code, reason) throw exc::domain_error(code, __FILE__, STRINGIZE(__LINE__), reason)
#define THROW_RUNTIME_ERROR(code, reason) throw exc::runtime_error(code, __FILE__, STRINGIZE(__LINE__), reason)
#define THROW_OVERFLOW_ERROR(code, reason) throw exc::overflow_error(code, __FILE__, STRINGIZE(__LINE__), reason)
#define THROW_RANGE_ERROR(code, reason) throw exc::range_error(code, __FILE__, STRINGIZE(__LINE__), reason)
#define THROW_UNDERFLOW_ERROR(code, reason) throw exc::underflow_error(code, __FILE__, STRINGIZE(__LINE__), reason)
