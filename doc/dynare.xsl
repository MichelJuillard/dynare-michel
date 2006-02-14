<?xml version='1.0'?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
                xmlns:fo="http://www.w3.org/1999/XSL/Format"
                xmlns:doc="http://nwalsh.com/xsl/documentation/1.0"
                exclude-result-prefixes="doc"
                version='1.0'>

<!-- It is important to use indent="no" here, otherwise verbatim -->
<!-- environments get broken by indented tags...at least when the -->
<!-- callout extension is used...at least with some processors -->

<xsl:param name="generate.index">1</xsl:param>
<xsl:param name="refentry.generate.name">0</xsl:param>
<xsl:param name="refentry.generate.title">1</xsl:param>
<xsl:param name="section.autolabel" select="1"></xsl:param>
<xsl:param name="biblioentry.item.separator">, </xsl:param>

<xsl:attribute-set name="section.level1.properties">
  <xsl:attribute name="break-before">page</xsl:attribute>
</xsl:attribute-set>


</xsl:stylesheet>
