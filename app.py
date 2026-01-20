"""
Entry point for HuggingFace Spaces - UAB Theme
"""

from sud_promise_uab_theme import create_interface

if __name__ == "__main__":
    demo = create_interface()
    demo.launch()

