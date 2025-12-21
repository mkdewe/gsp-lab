"""Logowanie"""

from config import DEFAULT_LOGGING

class Logger:
    _enabled = DEFAULT_LOGGING
    
    @classmethod
    def set_enabled(cls, enabled): 
        cls._enabled = enabled
    
    @staticmethod
    def info(msg): 
        if Logger._enabled: print(f"INFO: {msg}")
    
    @staticmethod
    def error(msg): 
        print(f"ERROR: {msg}")
    
    @staticmethod
    def section(title): 
        if Logger._enabled: print(f"\n=== {title} ===")
    
    @staticmethod
    def debug(msg): 
        if Logger._enabled: print(f"DEBUG: {msg}")