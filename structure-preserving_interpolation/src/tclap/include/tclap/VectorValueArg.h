/****************************************************************************** 
 * 
 *  file:  VectorValueArg.h
 * 
 * 
 *  See the file COPYING in the top directory of this distribution for
 *  more information.
 *  
 *  THE SOFTWARE IS PROVIDED _AS IS_, WITHOUT WARRANTY OF ANY KIND, EXPRESS 
 *  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
 *  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
 *  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 *  DEALINGS IN THE SOFTWARE.  
 *  
 *****************************************************************************/

#ifndef TCLAP_VECTORVALUE_ARGUMENT_H
#define TCLAP_VECTORVALUE_ARGUMENT_H

#include <string>
#include <vector>

#include <tclap/Arg.h>
#include <tclap/Constraint.h>

namespace TCLAP
{

/**
 * The basic labeled argument that parses a value.
 * This is a template class, which means the type T defines the type
 * that a given object will attempt to parse when the flag/name is matched
 * on the command line.  While there is nothing stopping you from creating
 * an unflagged ValueArg, it is unwise and would cause significant problems.
 * Instead use an UnlabeledValueArg.
 */
template<class T>
class VectorValueArg: public Arg
{
protected:

  /**
   * The value parsed from the command line.
   * Can be of any type, as long as the >> operator for the type
   * is defined.
   */
  T _value;

  std::vector<T> _valueVector;
  /**
   * Used to support the reset() method so that ValueArg can be
   * reset to their constructed value.
   */
  T _default;

  /**
   * A human readable description of the type to be parsed.
   * This is a hack, plain and simple.  Ideally we would use RTTI to
   * return the name of type T, but until there is some sort of
   * consistent support for human readable names, we are left to our
   * own devices.
   */
  std::string _typeDesc;

  /**
   * A Constraint this Arg must conform to.
   */
  Constraint<T>* _constraint;

  /**
   * Extracts the values from the string.
   * Attempts to parse string as type T or as vector<T>, if this fails an exception
   * is thrown.
   * \param val - value to be parsed.
   */
  void _extractVectorValue( std::vector<T>& valVec, const std::string& argStr );

  /**
   * Check wheteher the string contains the next argument flag or argument name
   */
  bool isArgumentFlag( const std::string& val ) const;

public:

  /**
   * Labeled ValueArg constructor.
   * You could conceivably call this constructor with a blank flag,
   * but that would make you a bad person.  It would also cause
   * an exception to be thrown.   If you want an unlabeled argument,
   * use the other constructor.
   * \param flag - The one character flag that identifies this
   * argument on the command line.
   * \param name - A one word name for the argument.  Can be
   * used as a long flag on the command line.
   * \param desc - A description of what the argument is for or
   * does.
   * \param req - Whether the argument is required on the command
   * line.
   * \param value - The default value assigned to this argument if it
   * is not present on the command line.
   * \param typeDesc - A short, human readable description of the
   * type that this object expects.  This is used in the generation
   * of the USAGE statement.  The goal is to be helpful to the end user
   * of the program.
   * \param v - An optional visitor.  You probably should not
   * use this unless you have a very good reason.
   */
  VectorValueArg( const std::string& flag, const std::string& name, const std::string& desc, bool req, T value, const std::string& typeDesc, Visitor* v = NULL );

  /**
   * Labeled ValueArg constructor.
   * You could conceivably call this constructor with a blank flag,
   * but that would make you a bad person.  It would also cause
   * an exception to be thrown.   If you want an unlabeled argument,
   * use the other constructor.
   * \param flag - The one character flag that identifies this
   * argument on the command line.
   * \param name - A one word name for the argument.  Can be
   * used as a long flag on the command line.
   * \param desc - A description of what the argument is for or
   * does.
   * \param req - Whether the argument is required on the command
   * line.
   * \param value - The default value assigned to this argument if it
   * is not present on the command line.
   * \param typeDesc - A short, human readable description of the
   * type that this object expects.  This is used in the generation
   * of the USAGE statement.  The goal is to be helpful to the end user
   * of the program.
   * \param parser - A CmdLine parser object to add this Arg to
   * \param v - An optional visitor.  You probably should not
   * use this unless you have a very good reason.
   */
  VectorValueArg( const std::string& flag, const std::string& name, const std::string& desc, bool req, T value, const std::string& typeDesc, CmdLineInterface& parser, Visitor* v = NULL );

  /**
   * Labeled ValueArg constructor.
   * You could conceivably call this constructor with a blank flag,
   * but that would make you a bad person.  It would also cause
   * an exception to be thrown.   If you want an unlabeled argument,
   * use the other constructor.
   * \param flag - The one character flag that identifies this
   * argument on the command line.
   * \param name - A one word name for the argument.  Can be
   * used as a long flag on the command line.
   * \param desc - A description of what the argument is for or
   * does.
   * \param req - Whether the argument is required on the command
   * line.
   * \param value - The default value assigned to this argument if it
   * is not present on the command line.
   * \param constraint - A pointer to a Constraint object used
   * to constrain this Arg.
   * \param parser - A CmdLine parser object to add this Arg to.
   * \param v - An optional visitor.  You probably should not
   * use this unless you have a very good reason.
   */
  VectorValueArg( const std::string& flag, const std::string& name, const std::string& desc, bool req, T value, Constraint<T>* constraint, CmdLineInterface& parser, Visitor* v = NULL );

  /**
   * Labeled ValueArg constructor.
   * You could conceivably call this constructor with a blank flag,
   * but that would make you a bad person.  It would also cause
   * an exception to be thrown.   If you want an unlabeled argument,
   * use the other constructor.
   * \param flag - The one character flag that identifies this
   * argument on the command line.
   * \param name - A one word name for the argument.  Can be
   * used as a long flag on the command line.
   * \param desc - A description of what the argument is for or
   * does.
   * \param req - Whether the argument is required on the command
   * line.
   * \param value - The default value assigned to this argument if it
   * is not present on the command line.
   * \param constraint - A pointer to a Constraint object used
   * to constrain this Arg.
   * \param v - An optional visitor.  You probably should not
   * use this unless you have a very good reason.
   */
  VectorValueArg( const std::string& flag, const std::string& name, const std::string& desc, bool req, T value, Constraint<T>* constraint, Visitor* v = NULL );

  /**
   * Handles the processing of the argument.
   * This re-implements the Arg version of this method to set the
   * _value of the argument appropriately.  It knows the difference
   * between labeled and unlabeled.
   * \param i - Pointer the the current argument in the list.
   * \param args - Mutable list of strings. Passed
   * in from main().
   */
  virtual bool processArg( int* i, std::vector<std::string>& args );

  /**
   * Returns the value of the argument.
   */
  std::vector<T>& getValue();

  /**
   * Specialization of shortID.
   * \param val - value to be used.
   */
  virtual std::string shortID( const std::string& val = "val" ) const;

  /**
   * Specialization of longID.
   * \param val - value to be used.
   */
  virtual std::string longID( const std::string& val = "val" ) const;

  virtual void reset();

private:
  /**
   * Prevent accidental copying
   */
  VectorValueArg<T>( const VectorValueArg<T>& rhs );
  VectorValueArg<T>& operator=( const VectorValueArg<T>& rhs );
};

/**
 * Constructor implementation.
 */
template<class T>
VectorValueArg<T>::VectorValueArg( const std::string& flag, const std::string& name, const std::string& desc, bool req, T val, const std::string& typeDesc, Visitor* v ) :
    Arg( flag, name, desc, req, true, v ), _value( val ), _valueVector( 1, val ), _default( val ), _typeDesc( typeDesc ), _constraint( NULL )
{
}

template<class T>
VectorValueArg<T>::VectorValueArg( const std::string& flag, const std::string& name, const std::string& desc, bool req, T val, const std::string& typeDesc, CmdLineInterface& parser, Visitor* v ) :
    Arg( flag, name, desc, req, true, v ), _value( val ), _valueVector( 1, val ), _default( val ), _typeDesc( typeDesc ), _constraint( NULL )
{
  parser.add( this );
}

template<class T>
VectorValueArg<T>::VectorValueArg( const std::string& flag, const std::string& name, const std::string& desc, bool req, T val, Constraint<T>* constraint, Visitor* v ) :
    Arg( flag, name, desc, req, true, v ), _value( val ), _valueVector( 1, val ), _default( val ), _typeDesc( constraint->shortID() ), _constraint( constraint )
{
}

template<class T>
VectorValueArg<T>::VectorValueArg( const std::string& flag, const std::string& name, const std::string& desc, bool req, T val, Constraint<T>* constraint, CmdLineInterface& parser, Visitor* v ) :
    Arg( flag, name, desc, req, true, v ), _value( val ), _valueVector( 1, val ), _default( val ), _typeDesc( constraint->shortID() ), _constraint( constraint )
{
  parser.add( this );
}

/**
 * Implementation of getValue().
 */
template<class T>
std::vector<T>& VectorValueArg<T>::getValue()
{
  return _valueVector;
}

/**
 * Implementation of processArg().
 */
template<class T>
bool VectorValueArg<T>::processArg( int *i, std::vector<std::string>& args )
{
  if( _ignoreable && Arg::ignoreRest() )
    return false;

  if( _hasBlanks( args[*i] ) )
    return false;

  std::string flag = args[*i];

  std::string value = "";
  trimFlag( flag, value );

  if( argMatches( flag ) )
  {
    if( _alreadySet )
    {
      if( _xorSet )
        throw(CmdLineParseException( "Mutually exclusive argument already set!", toString() ));
      else
        throw(CmdLineParseException( "Argument already set!", toString() ));
    }

    if( Arg::delimiter() != ' ' && value == "" )
      throw(ArgParseException( "Couldn't find delimiter for this argument!", toString() ));

    // clear default values from vector
    _valueVector.clear();

    if( value == "" )
    {
        while( static_cast<unsigned int>( *i ) < (args.size()-1) )
        {
          if( isArgumentFlag( args[*i+1] ) ) break;
          (*i)++;

          std::string currentArg = args[*i];
          std::vector<T> currentValueVector;

          _extractVectorValue( currentValueVector, currentArg );

          // append parsed values
          _valueVector.insert(_valueVector.end(), currentValueVector.begin(), currentValueVector.end());
        }
        if( _valueVector.size() == 0 )
                throw(ArgParseException( "Missing a value(s) for this argument!", toString() ));
     }
    else
      _extractVectorValue( _valueVector, value );

    _alreadySet = true;
    _checkWithVisitor();
    return true;
  }
  else
    return false;
}

/**
 * Implementation of shortID.
 */
template<class T>
std::string VectorValueArg<T>::shortID( const std::string& val ) const
{
  static_cast<void>( val ); // Ignore input, don't warn
  return Arg::shortID( _typeDesc );
}

/**
 * Implementation of longID.
 */
template<class T>
std::string VectorValueArg<T>::longID( const std::string& val ) const
{
  static_cast<void>( val ); // Ignore input, don't warn
  return Arg::longID( _typeDesc );
}

template<class T>
bool VectorValueArg<T>::isArgumentFlag( const std::string& val ) const
{
  if( val.length() < 2 )
    return false;

  if( val.substr( 0, Arg::flagStartString().length() ) == Arg::flagStartString() || val.substr( 0, Arg::nameStartString().length() ) == Arg::nameStartString() )
  {
    // check for negative number
    if( "-" == Arg::flagStartString() && val[0] == '-' )
    {
      std::string number = "0123456789";
      static const std::basic_string<char>::size_type npos = std::basic_string<char>::npos;
      if( npos != number.find( val.substr( 1, 1 ) ) )
        return false;
    }

    return true;
  }

  return false;
}

/*
 * Extract a value of type T from it's string representation contained
 * in strVal. The ValueLike parameter used to select the correct
 * specialization of ExtractValue depending on the value traits of T.
 * StringLike uses assignment (operator=) to assign from strVal.
 */
template<class T>
void SetVectorFromString(std::vector<T> &v, const std::string &s)
{
    std::istringstream iss(s);
    while (iss.good()) {
      if ( iss.peek() != EOF )
      {
        T tmp;
        iss >> tmp;
        v.push_back(tmp);
      }
    }
    if ( iss.fail() )
      throw( ArgParseException("Couldn't read argument value "
         "from string '" + s + "'"));
}

template<class T>
void VectorValueArg<T>::_extractVectorValue( std::vector<T>& valVec, const std::string& argStr )
{
    try
    {
      SetVectorFromString(valVec, argStr);

      //    std::cout<<"\n val="<<val<<" -- "<<(val.find(' ') != std::string::npos)<<"\n";
  //    if(val.find(' ') != std::string::npos)
  //      ExtractValue( _valueVector, val, StringLike() );
  //    else
  //      ExtractValue( _value, val, typename ArgTraits<T>::ValueCategory() );
    }
    catch( ArgParseException &e )
    {
      throw ArgParseException( e.error(), toString() );
    }

    if( _constraint != NULL )
    {
      for ( int i = 0; static_cast<unsigned int>(i) < valVec.size(); i++ )
        if( !_constraint->check( valVec[i] ) )
          throw(CmdLineParseException( "Value '" + argStr +"' does not meet constraint: " + _constraint->description(), toString() ));
    }

}

template<class T>
void VectorValueArg<T>::reset()
{
  Arg::reset();
  _value = _default;
  _valueVector.clear();
  _valueVector.push_back(_value);
}

} // namespace TCLAP

#endif
