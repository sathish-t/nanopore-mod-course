---
layout: page
title: Exercises
---


<table>
  <tr>
    <th>Topic</th>
    <th>Title</th>
    <th>Output</th>
  </tr>
{% for exercise in site.pages %}
  {% if exercise.layout == 'exercise' %}
    {% unless exercise.title contains "Template" %}
      <tr>
      <td nowrap>{{ exercise.topic | replace:'and','&'  }}</td>
      <td nowrap><a href="{{ exercise.url | prepend: site.baseurl }}">
        {{ exercise.title }}</a></td>
      {% capture output_file %}{{ exercise.url | remove: 'exercises' | remove: '/' | prepend: '/solutions/' }}{% endcapture %}
      <td>
      {% for solution in site.static_files %}
        {% if solution.path contains output_file %}
          <a href="{{ solution.path | prepend: site.baseurl}}">
            [{{ solution.path | split:"." | last}}]</a>
        {% endif %}
      {% endfor %}
      </td>
      </tr>
    {% endunless %}
  {% endif %}
{% endfor %}
</table>
<p align="right"><font size="-1">
  <a href="{{ site.baseurl }}/exercises/">[back to top]</a>
</font></p>