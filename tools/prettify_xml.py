import xml.etree.ElementTree as eTree
import xml.dom.minidom as minidom

# Read the XML file
tree        =   eTree.parse('../src/JobFile.xml')
root        =   tree.getroot()
xmlSting    =   eTree.tostring(root, encoding = "utf-8", method = "xml")
dom         =   minidom.parseString(xmlSting)
prettyXML   =   dom.toprettyxml()

with open(f"test.xml", "w") as outObject:
    outObject.write(prettyXML)