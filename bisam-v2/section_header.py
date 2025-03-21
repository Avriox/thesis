#!/usr/bin/env python3
import sys
import os
import re


def process_file(file_path):
    """Process the file and replace both section and subsection markers"""
    with open(file_path, 'r') as f:
        lines = f.readlines()

    modified_lines = []
    modified = False
    line_number = 0

    for line in lines:
        line_number += 1
        if "//section" in line:
            # Extract section name
            match = re.search(r'^\s*//section\s+(.+?)(?:\s*$)', line)
            if match:
                section = match.group(1).strip()
                modified = True
                modified_lines.extend(generate_section_header(section, file_path, "h1"))
                print(f"Line {line_number}: Replaced section marker with header for '{section}'")
            else:
                modified_lines.append(line)
        elif "//subsection" in line:
            # Extract subsection name
            match = re.search(r'^\s*//subsection\s+(.+?)(?:\s*$)', line)
            if match:
                section = match.group(1).strip()
                modified = True
                modified_lines.extend(generate_section_header(section, file_path, "h2"))
                print(f"Line {line_number}: Replaced subsection marker with header for '{section}'")
            else:
                modified_lines.append(line)
        else:
            modified_lines.append(line)

    if modified:
        with open(file_path, 'w') as f:
            f.writelines(modified_lines)
        print(f"File modified: {file_path}")
        return True

    print(f"No section markers found in {file_path}")
    return False


def generate_section_header(section, file_path, level):
    """Generate a formatted section header"""
    # Detect language from file extension
    file_ext = os.path.splitext(file_path)[1].lower()
    comment_config = {
        ".cpp": ("/*", "*/", 80),  # C/C++
        ".c": ("/*", "*/", 80),  # C
        ".h": ("/*", "*/", 80),  # C/C++ header
        ".hpp": ("/*", "*/", 80),  # C++ header
        ".py": ("#", "", 80),  # Python
        ".html": ("<!--", "-->", 80),  # HTML
        ".js": ("/*", "*/", 80),  # JavaScript
        ".css": ("/*", "*/", 80),  # CSS
        ".java": ("/*", "*/", 80),  # Java
        ".kt": ("/*", "*/", 80),  # Kotlin
        ".ts": ("/*", "*/", 80),  # TypeScript
    }

    comment_start, comment_end, line_length = comment_config.get(file_ext, ("//", "", 80))

    # Calculate available space for dashes (adjust for comment syntax and text)
    if level == "h1":
        # Main header (3 lines)
        available_length = line_length - len(comment_start) - len(comment_end) - 2  # Account for spaces
        dashes = "-" * available_length

        # For second line with text, center text with spaces
        text_length = len(section)
        available_for_spaces = available_length - text_length
        spaces_left = " " * (available_for_spaces // 2)
        spaces_right = " " * (available_for_spaces // 2 + (available_for_spaces % 2))

        output = [
            f"{comment_start}{' ' if comment_start[-1] != ' ' else ''}{dashes}{' ' if comment_end and comment_end[0] != ' ' else ''}{comment_end}\n",
            f"{comment_start}{' ' if comment_start[-1] != ' ' else ''}{spaces_left}{section}{spaces_right}{' ' if comment_end and comment_end[0] != ' ' else ''}{comment_end}\n",
            f"{comment_start}{' ' if comment_start[-1] != ' ' else ''}{dashes}{' ' if comment_end and comment_end[0] != ' ' else ''}{comment_end}\n"
        ]
    elif level == "h2":
        # Subheader (single line)
        available_length = line_length - len(comment_start) - len(comment_end) - len(section) - 4  # Just 2 spaces
        dashes_left = "-" * (available_length // 2)
        dashes_right = "-" * (available_length // 2 + (available_length % 2))
        output = [
            f"{comment_start}{' ' if comment_start[-1] != ' ' else ''}{dashes_left} {section} {dashes_right}{' ' if comment_end and comment_end[0] != ' ' else ''}{comment_end}\n"]
    else:
        return []  # Return empty list if level is invalid

    return output


def main():
    # Check arguments
    if len(sys.argv) < 2:
        print("Usage: section_header.py <file_path>")
        sys.exit(1)

    file_path = sys.argv[1]

    if not os.path.exists(file_path):
        print(f"Error: File {file_path} does not exist")
        sys.exit(1)

    result = process_file(file_path)
    if result:
        print(f"Successfully processed {file_path}")
    else:
        print(f"No changes made to {file_path}")


if __name__ == "__main__":
    main()
