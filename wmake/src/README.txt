Minimal changes to lemon for C++ integration.

- New '-e' command line option to define the code extension.
  By default this is 'c', but with this option can define something
  like -ecxx etc for using a C++ compiler.

- New '%namespace' directive. This can be used to embed the 'Parse*'
  routines into a C++ namespace. The namespace can be anonymous or
  contain multiple nested namespaces. For example,

      %namespace {}
      %namespace {ns1::ns2::ns3}


  One simple means to encapsulate code is to use an anonymous
  namespace for the Lemon code and place all C++ interface code with
  the %code block in the same translation unit.

  This allows good encapsulation without fundamentally changing how
  Lemon works.

--
2019-09-27
