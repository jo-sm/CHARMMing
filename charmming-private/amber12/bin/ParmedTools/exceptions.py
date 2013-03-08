""" Exceptions used in Parmed """

class ParmError(Exception):
   """ Base parmed error """
   def __init__(self, msg='parmed error'):
      self.msg = msg
   def __str__(self):
      return self.msg

class ChangeRadiiError(ParmError):
   pass

class WriteOFFError(ParmError):
   pass

class ParmedUtilsError(ParmError):
   pass

class ParmedChangeError(ParmError):
   pass

class ParmedAddLJTypeError(ParmError):
   pass

class ChangeLJPairError(ParmError):
   pass

class LJ_TypeError(ParmError):
   pass

class ParmedMoleculeError(ParmError):
   pass

class CoarseGrainError(ParmError):
   pass

class ChangeStateError(ParmError):
   pass

class SetParamError(ParmError):
   pass

class DeleteDihedralError(ParmError):
   pass
