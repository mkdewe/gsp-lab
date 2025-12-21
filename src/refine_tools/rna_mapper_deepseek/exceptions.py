"""Wyjątki"""

class RNAMapperError(Exception): pass
class StructureParseError(RNAMapperError): pass
class PDBFixError(RNAMapperError): pass